# trialforge monorepo: the Rust core lives in crates/, the Python package in bindings/python and
# the R package in bindings/r. Both bindings replay the shared conformance cases in spec/.
# Targets name their language (`rs-*`, `py-*`, `r-*`), and the bare names (`test`, `check`,
# `install`) run all three.

PYTHON ?= $(CURDIR)/.venv/bin/python
RSCRIPT ?= Rscript
PY_DIR := bindings/python
R_DIR := bindings/r

.PHONY: help
help: ## Show this help message
	@echo 'Usage: make [target]'
	@echo ''
	@echo 'Available targets:'
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "  %-20s %s\n", $$1, $$2}' $(MAKEFILE_LIST)

# -- All ----------------------------------------------------------------------------------------

.PHONY: install
install: py-install r-install ## Build and install the Python and R packages for development

.PHONY: test
test: rs-test py-test r-test ## Run the Rust, Python and R test suites

.PHONY: check
check: rs-check py-check r-test ## Run every lint, type check and test suite (the pre-push gate)

# -- Rust ---------------------------------------------------------------------------------------

.PHONY: rs-test
rs-test: ## Run the Rust workspace tests
	cargo test --workspace --exclude trialforge-python --all-features

.PHONY: rs-lint
rs-lint: ## Format the Rust code and apply clippy fixes
	cargo fmt --all
	cargo clippy --workspace --exclude trialforge-python --all-targets --all-features --fix --allow-dirty

.PHONY: rs-check
rs-check: ## Rust gate: rustfmt check, clippy (warnings are errors), tests
	cargo fmt --all --check
	cargo fmt --manifest-path $(R_DIR)/src/rust/Cargo.toml --check
	cargo clippy --workspace --exclude trialforge-python --all-targets --all-features -- -D warnings
	cargo clippy -p trialforge-python -- -D warnings
	cargo clippy --manifest-path $(R_DIR)/src/rust/Cargo.toml --lib -- -D warnings
	cargo test --workspace --exclude trialforge-python --all-features

# -- Python -------------------------------------------------------------------------------------

.venv:
	uv venv --python 3.12 .venv

.PHONY: py-install
py-install: .venv ## Build the extension and install the Python package (editable) into .venv
	uv pip install --python "$(PYTHON)" "maturin>=1.8,<2" pytest ruff pyright polars pandas pyarrow pyreadstat "pointblank>=1.1.1" greenwood lxml xmlschema odmlib
	cd $(PY_DIR) && VIRTUAL_ENV=$(CURDIR)/.venv $(CURDIR)/.venv/bin/maturin develop --uv

.PHONY: py-test
py-test: ## Run the Python test suite (run py-install first after Rust changes)
	cd $(PY_DIR) && $(PYTHON) -m pytest

.PHONY: py-lint
py-lint: ## Run the ruff formatter and linter (with fixes)
	cd $(PY_DIR) && $(PYTHON) -m ruff format
	cd $(PY_DIR) && $(PYTHON) -m ruff check --fix

.PHONY: py-type-check
py-type-check: ## Type-check the Python package (pyright strict)
	cd $(PY_DIR) && $(PYTHON) -m pyright

.PHONY: py-check
py-check: ## Python gate: ruff, pyright, pytest
	cd $(PY_DIR) && $(PYTHON) -m ruff format --check
	cd $(PY_DIR) && $(PYTHON) -m ruff check
	cd $(PY_DIR) && $(PYTHON) -m pyright
	cd $(PY_DIR) && $(PYTHON) -m pytest

# -- R ------------------------------------------------------------------------------------------

# `R CMD build` copies the package away from the repo, which breaks the path dependency on
# crates/trialforge, so r-build first stages a self-contained copy with the engine crates and
# vendored third-party crates (scripts/stage_r_package.py, into build/r-package/).

.PHONY: r-document
r-document: ## Compile the Rust library and regenerate R wrappers, NAMESPACE and man pages
	@# pkgbuild only watches src/, not the core crates: drop the objects so it always relinks.
	rm -f $(R_DIR)/src/*.so $(R_DIR)/src/*.o
	$(RSCRIPT) -e 'devtools::document("$(R_DIR)")'

.PHONY: r-test
r-test: ## Run the R test suite
	$(RSCRIPT) -e 'devtools::test("$(R_DIR)", stop_on_failure = TRUE)'

.PHONY: r-build
r-build: ## Stage a self-contained R package (engine crates + vendored deps) and build the tarball
	$(PYTHON) scripts/stage_r_package.py
	cd build/r-package && rm -f trialforge_*.tar.gz && R CMD build --no-manual trialforge

.PHONY: r-check
r-check: r-build ## R CMD check --as-cran on the staged tarball (offline build from vendored crates)
	cd build/r-package && _R_CHECK_FORCE_SUGGESTS_=false R CMD check --as-cran --no-manual trialforge_*.tar.gz

# admiral equivalence: admiral re-derives our ADaM from our SDTM; every compared value must match.
# Runs in an isolated `ir` environment (nothing is installed into the user library).
ADMIRAL_IR := ir run --isolated --r-version 4.6 --exclude-newer 2026-09-25 --with admiral --with admiralonco --with admiralophtha --with admiralpeds --with admiralvaccine --with haven --with dplyr --with stringr --with tidyr
PKNCA_IR := ir run --isolated --r-version 4.6 --exclude-newer 2026-09-25 --with PKNCA --with haven --with dplyr
ADMIRAL_CASES := parallel_continuous:{} \
	'parallel_continuous:{"indication":"hypertension","weeks":12,"dropout_rate":0.9}' \
	'parallel_continuous:{"indication":"obesity","ae_rate_multiplier":4}' \
	oncology_solid_tumor:{} \
	'oncology_solid_tumor:{"orr":[0.6,0.3],"cutoff_months":22,"withdrawal_prob":0.3,"os_pfs_correlation":0.9}' \
	phase1_sad_mad:{} \
	'phase1_sad_mad:{"part":"mad","doses":[50,100,200,400],"compartments":2}' \
	crossover_be:{} \
	'crossover_be:{"design":"2x2x4","dropout_rate":3}' \
	'crossover_be:{"design":"food_effect","gmr":1.3,"ka_ratio":0.5}' \
	dose_ranging_ph2:{} \
	'parallel_continuous:{"comparator":"active","ni_margin":0.3,"indication":"obesity"}' \
	cns_alzheimers:{} \
	'cns_alzheimers:{"item_missing":0.15,"dropout_rate":2,"weeks":48}' \
	cv_outcomes:{} \
	'cv_outcomes:{"target_events":40,"withdrawal_rate":0.3,"discontinuation_rate":0.3,"components":[1,1,1]}' \
	ophthalmology:{} \
	'ophthalmology:{"interval_weeks":[4,16],"loading_doses":1,"weeks":24,"dropout_rate":2,"missed_visit":0.3}' \
	pediatric_growth:{} \
	'pediatric_growth:{"age_range":[3,15],"weeks":104,"dropout_rate":1,"height_velocity":[6,12]}' \
	respiratory_exacerbation:{} \
	'respiratory_exacerbation:{"exacerbation_rate":3,"dispersion":2,"severe_fraction":0.6,"weeks":24,"dropout_rate":2}' \
	vaccine:{} \
	'vaccine:{"doses":1,"reactogenicity":2.5,"fever_rate":0.5,"seropositive":0.5,"dropout_rate":2,"missed_visit":0.3}' \
	'oncology_solid_tumor:{"crossover_prob":0.7}' \
	'oncology_solid_tumor:{"interim_months":12}' \
	'cv_outcomes:{"target_events":40,"interim_events":20}' \
	heme_lymphoma:{} \
	'heme_lymphoma:{"orr":[0.8,0.2],"cutoff_months":22,"withdrawal_prob":0.3,"os_pfs_correlation":0.9}' \
	'heme_lymphoma:{"interim_months":9}' \
	dose_escalation:{} \
	'dose_escalation:{"method":"3+3","dropout_rate":2}' \
	'dose_escalation:{"method":"mtpi2","cohort_size":2,"weeks":4}'

.PHONY: r-admiral
r-admiral: ## Compare our ADaM with admiral's derivations from our SDTM (needs `ir`)
	@for case in $(ADMIRAL_CASES); do \
		archetype=$${case%%:*}; controls=$${case#*:}; out=build/admiral/$$archetype; \
		rm -rf $$out && $(PYTHON) -c "import json, sys, trialforge as tf; tf.study(sys.argv[1], n=300, seed=21, **json.loads(sys.argv[2])).write_submission(sys.argv[3])" "$$archetype" "$$controls" $$out && \
		echo "== $$archetype $$controls" && \
		$(ADMIRAL_IR) scripts/admiral_equivalence.R $$out 2>&1 | grep -E "admiral [0-9]|differs|rows" || exit 1; \
		if [ -f $$out/sdtm/pp.xpt ]; then $(PKNCA_IR) scripts/pknca_equivalence.R $$out 2>&1 | grep -E "PKNCA [0-9]|differ" || exit 1; fi; \
	done

.PHONY: r-install
r-install: ## Install the R package from the source tree (in place, so the core path resolves)
	R CMD INSTALL --no-multiarch $(R_DIR)

# -- Housekeeping -------------------------------------------------------------------------------

.PHONY: clean
clean: ## Remove Rust, Python and R build artifacts
	cargo clean
	rm -rf $(R_DIR)/src/rust/target $(R_DIR)/src/*.o $(R_DIR)/src/*.so $(R_DIR)/src/*.dll
	rm -f $(PY_DIR)/python/trialforge/_core*.so $(PY_DIR)/python/trialforge/_core*.pyd

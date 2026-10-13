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
	cargo test --workspace --exclude trialforge-python

.PHONY: rs-lint
rs-lint: ## Format the Rust code and apply clippy fixes
	cargo fmt --all
	cargo clippy --workspace --exclude trialforge-python --all-targets --fix --allow-dirty

.PHONY: rs-check
rs-check: ## Rust gate: rustfmt check, clippy (warnings are errors), tests
	cargo fmt --all --check
	cargo fmt --manifest-path $(R_DIR)/src/rust/Cargo.toml --check
	cargo clippy --workspace --exclude trialforge-python --all-targets -- -D warnings
	cargo test --workspace --exclude trialforge-python

# -- Python -------------------------------------------------------------------------------------

.venv:
	uv venv --python 3.12 .venv

.PHONY: py-install
py-install: .venv ## Build the extension and install the Python package (editable) into .venv
	uv pip install --python "$(PYTHON)" "maturin>=1.8,<2" pytest ruff pyright polars pandas pyarrow pyreadstat
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

# There is no r-check target yet. `R CMD build` copies the package away from the repo, which
# breaks the path dependency on crates/trialforge. Vendoring the core crates fixes that.

.PHONY: r-document
r-document: ## Compile the Rust library and regenerate R wrappers, NAMESPACE and man pages
	$(RSCRIPT) -e 'devtools::document("$(R_DIR)")'

.PHONY: r-test
r-test: ## Run the R test suite
	$(RSCRIPT) -e 'devtools::test("$(R_DIR)", stop_on_failure = TRUE)'

.PHONY: r-install
r-install: ## Install the R package from the source tree (in place, so the core path resolves)
	R CMD INSTALL --no-multiarch $(R_DIR)

# -- Housekeeping -------------------------------------------------------------------------------

.PHONY: clean
clean: ## Remove Rust, Python and R build artifacts
	cargo clean
	rm -rf $(R_DIR)/src/rust/target $(R_DIR)/src/*.o $(R_DIR)/src/*.so $(R_DIR)/src/*.dll
	rm -f $(PY_DIR)/python/trialforge/_core*.so $(PY_DIR)/python/trialforge/_core*.pyd

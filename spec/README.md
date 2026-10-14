# spec/

Language-neutral contracts shared by the Rust core and both bindings. If behaviour is defined
here, the Python and R packages must agree on it exactly.

- `study-spec.schema.json`: the study spec. This is v0, a draft.
- `conformance/`: golden cases that both test suites replay.
  - `rng_streams.json`: RNG stream seeds (`crates/trialforge-core/src/rng.rs`). Seeds are
    decimal strings because R has no unsigned 64-bit integer.
  - `study_xpt.json`: for reference `study(archetype, n, seed, **controls)` calls, the row count of each
    dataset, the MD5 of each XPT file written and of each define.xml. Both bindings write through the Rust writer
    and hash the bytes themselves, which proves Python and R produce identical files.
- `designs/`: per-archetype semantics, i.e. controls, defaults, derived truths and expected
  conformance warnings. Empty so far.
- `derivations/`: ADaM derivation rules, each with the admiral function it mirrors. These rules
  are also the source of the define.xml Methods. Empty so far.

Changing a golden value changes generated data. Do it on purpose, update every affected case,
and say so in the commit message.

To regenerate `study_xpt.json` after a deliberate change (run from the repo root after
`make py-install`):

```bash
.venv/bin/python - <<'PY'
import hashlib, json, pathlib, tempfile
import trialforge as tf
path = pathlib.Path("spec/conformance/study_xpt.json")
doc = json.loads(path.read_text())
for case in doc["cases"]:
    s = tf.study(case["archetype"], n=case["n"], seed=case["seed"], **case["controls"])
    case["rows"] = {k: len(v) for k, v in {**s.sdtm, **s.adam}.items()}
    out = pathlib.Path(tempfile.mkdtemp())
    case["xpt_md5"] = {
        p.name: hashlib.md5(p.read_bytes()).hexdigest() for p in s.write(out) if p.suffix == ".xpt"
    }
    case["define_md5"] = {
        std: hashlib.md5(s.define_xml(std).encode()).hexdigest() for std in ("sdtm", "adam")
    }
path.write_text(json.dumps(doc, indent=2) + "\n")
PY
```

## Language adaptations

Package-level names match across Python and R. Where they cannot, the difference is recorded
here.

| Python | R | Why |
|---|---|---|
| `Study.write()` | `write_study()` | `write()` would mask `base::write()` |
| `Study.verify(level=3)` | not available (`verify(trial, level = 3)` errors) | Level 3 runs Pointblank's CDISC conformance engine, which only exists in the Python package; the R `pointblank` package has no SDTMIG rules |
| `Study.write(formats=["parquet"])` | not available | The Parquet writer is compiled only into the Python build, to keep the R package's vendored crates small |
| Date-times (e.g. ADPC `ADTM`) are naive timestamps (no time zone) | `POSIXct` with `tzone = "UTC"` | The study clock has no time zone; R needs one, and UTC keeps the clock time unchanged |


# trialforge

Synthetic CDISC clinical-trial data from declarative study specs: trial design datasets, SDTM,
ADaM and define.xml for a complete mock study, with interpretable controls for many study
types and output that is verified to be workable.

trialforge has a Rust core with Python and R bindings. The same spec and seed produce identical
data in both languages.

> **Status: pre-alpha.** Nothing generates a study yet.

## Layout

| Path | What |
|---|---|
| `crates/` | Rust workspace. `trialforge` is the facade crate; `trialforge-*` are the engine crates. |
| `bindings/python/` | Python package (pyo3 + maturin). |
| `bindings/r/` | R package (extendr). |
| `spec/` | Language-neutral contracts: the study-spec JSON Schema and golden conformance cases that both bindings replay. |

## Development

Requirements: Rust ≥ 1.85, [uv](https://docs.astral.sh/uv/), and R ≥ 4.2 with `devtools` and
`rextendr`.

```bash
make py-install   # build the extension into .venv (re-run after Rust changes)
make r-install    # build and install the R package from the source tree
make check        # Rust, Python and R gates
make help         # all targets
```

## License

MIT

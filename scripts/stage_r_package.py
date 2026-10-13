"""Stage a self-contained copy of the R package for `R CMD build` / `R CMD check`.

The R package's Rust crate depends on `crates/trialforge` by relative path, which only resolves
inside this repository. This script builds `build/r-package/trialforge/` where:

- `src/rust/crates/` holds copies of the engine crates the R crate needs (not the CLI or the
  Python binding),
- `src/rust/Cargo.toml` becomes the root of a small Cargo workspace over those crates, with the
  repository's shared `[workspace.package]`, `[workspace.dependencies]` and lints,
- the Arrow/Parquet dependencies are stripped from `trialforge-io` (the R build never enables
  them, so they need not be vendored),
- third-party crates are vendored into `src/rust/vendor.tar.xz` with `vendor-config.toml`, so
  the package builds offline as CRAN requires.

Run from the repository root (needs network access, or a warm cargo cache, for `cargo vendor`):

    python scripts/stage_r_package.py
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STAGE = ROOT / "build" / "r-package" / "trialforge"
SKIP_CRATES = {"trialforge-cli"}


def section(text: str, header: str) -> str:
    """The body of a `[header]` TOML table (up to the next table)."""
    m = re.search(rf"^\[{re.escape(header)}\]\n(.*?)(?=^\[|\Z)", text, re.S | re.M)
    if not m:
        raise SystemExit(f"[{header}] not found in the workspace Cargo.toml")
    return m.group(1).strip()


def main() -> int:
    if STAGE.exists():
        shutil.rmtree(STAGE)
    ignore = shutil.ignore_patterns("target", "*.o", "*.so", "*.dll", ".cargo", "vendor", "Makevars")
    shutil.copytree(ROOT / "bindings" / "r", STAGE, ignore=ignore)
    rust = STAGE / "src" / "rust"
    for crate in sorted((ROOT / "crates").iterdir()):
        if crate.is_dir() and crate.name not in SKIP_CRATES:
            shutil.copytree(crate, rust / "crates" / crate.name, ignore=shutil.ignore_patterns("target"))
    shutil.copy(ROOT / "LICENSE", STAGE / "LICENSE.note")

    # The R build never compiles Rust tests or examples: drop them and their dev-dependencies
    # (proptest alone would pull in ~40 MB of platform crates).
    for crate in (rust / "crates").iterdir():
        for extra in ("tests", "examples"):
            shutil.rmtree(crate / extra, ignore_errors=True)
        manifest = crate / "Cargo.toml"
        manifest.write_text(
            re.sub(r"^\[dev-dependencies\]\n(.*?)(?=^\[|\Z)", "", manifest.read_text(), flags=re.S | re.M)
        )

    # trialforge-io without the Arrow/Parquet dependencies (features kept, but empty).
    io = rust / "crates" / "trialforge-io" / "Cargo.toml"
    text = io.read_text()
    text = re.sub(r"^(arrow-array|arrow-schema|parquet) = \{ version.*\n", "", text, flags=re.M)
    text = re.sub(r'^arrow = \[.*\]$', "arrow = []", text, flags=re.M)
    text = re.sub(r'^parquet = \[.*\]$', 'parquet = ["arrow"]', text, flags=re.M)
    io.write_text(text)

    # The R crate becomes a workspace root over the copied crates.
    ws = (ROOT / "Cargo.toml").read_text()
    manifest = rust / "Cargo.toml"
    text = manifest.read_text()
    text = re.sub(
        r'^trialforge = \{ path = "[^"]+" \}$',
        'trialforge = { path = "crates/trialforge" }',
        text,
        flags=re.M,
    )
    lints = "\n\n".join(
        f"[{h}]\n{section(ws, h)}" for h in ("workspace.lints.rust", "workspace.lints.clippy")
    )
    text += (
        '\n[workspace]\nresolver = "3"\nmembers = ["crates/*"]\n\n'
        f"[workspace.package]\n{section(ws, 'workspace.package')}\n\n"
        f"[workspace.dependencies]\n{section(ws, 'workspace.dependencies')}\n\n{lints}\n"
    )
    manifest.write_text(text)

    # Vendor third-party crates for offline builds.
    vendor = rust / "vendor"
    subprocess.run(["cargo", "vendor", "--quiet", str(vendor)], cwd=rust, check=True, stdout=subprocess.DEVNULL)
    subprocess.run(["tar", "-cJf", "vendor.tar.xz", "vendor"], cwd=rust, check=True)
    shutil.rmtree(vendor)
    (rust / "vendor-config.toml").write_text(
        '[source.crates-io]\nreplace-with = "vendored-sources"\n\n'
        '[source.vendored-sources]\ndirectory = "vendor"\n'
    )
    size = (rust / "vendor.tar.xz").stat().st_size / 1e6
    print(f"staged {STAGE.relative_to(ROOT)} (vendor.tar.xz: {size:.1f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

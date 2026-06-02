# Development

`seqfold`'s folding/Tm engine is written in **Rust** and exposed to Python via
[PyO3](https://pyo3.rs) / [maturin](https://www.maturin.rs). Pure-Rust code
lives in `src/core/`; the PyO3 bindings (the `seqfold._core` extension module)
are in `src/python.rs`; the Python package (the `Struct` re-exports and the
argparse CLI) is in `python/seqfold/`.

## Prerequisites

- A [Rust toolchain](https://rustup.rs)
- Python 3.8+
- [maturin](https://www.maturin.rs): `pip install maturin`

## Build & install locally

```bash
python -m venv .venv && . .venv/bin/activate
maturin develop --release        # compile the Rust extension into the venv
```

> Note: `maturin develop` drops the compiled module at
> `python/seqfold/_core.*.so`. If you ever switch ABI settings, delete stale
> `python/seqfold/_core*.so` first — Python prefers a version-specific `.so`
> over the `abi3` one, which can silently load an old build.

## Test

```bash
make test          # cargo test + the Python unittest suite
# or individually:
cargo test --lib                                   # Rust unit tests (engine internals)
python -m unittest discover tests -p '*_test.py'   # Python public API + CLI
```

`tests/fold_examples.py` regenerates `examples/dna.csv` / `examples/rna.csv`
(seqfold vs UNAFold). `examples/known_structures.fasta` holds real structured
sequences used by the benchmark.

`scripts/smoke_test.py` checks the public API + CLI of *any* install against
known-good values — handy for verifying a wheel from CI or a published release,
not just a local build:

```bash
python scripts/smoke_test.py        # exits non-zero on any mismatch

# against a wheel from a CI run, in a clean env:
gh run download <run-id> -n wheels-macos-aarch64 -D /tmp/whl   # pick your platform
python -m venv /tmp/v && /tmp/v/bin/pip install /tmp/whl/*.whl
/tmp/v/bin/python scripts/smoke_test.py

# against PyPI after a release:
python -m venv /tmp/v && /tmp/v/bin/pip install seqfold
/tmp/v/bin/python scripts/smoke_test.py
```

## Benchmark

```bash
cargo run --release --example bench                       # parallel
RAYON_NUM_THREADS=1 cargo run --release --example bench   # single-threaded
```

## Regenerating the energy tables

`src/core/data.rs` is generated from the original Python energy parameters by
`codegen/gen_data.py`. It's committed, so you only need to regenerate it if the
parameters change (the original `seqfold/dna.py` / `rna.py` are in git history).

## Releasing

**The GitHub Release tag is the single source of truth for the version.** You do
not edit a version file. The number in `Cargo.toml` is only a dev default; CI
overwrites it from the release tag before building.

To cut a release:

1. Make sure `main` is green.
2. Create a GitHub Release with the version as its tag:
   ```bash
   gh release create 0.11.0 --generate-notes        # use your tag convention
   ```
   (or via the GitHub UI: Releases → Draft a new release → create the tag.)

That's it. Publishing the release triggers `.github/workflows/release.yml`,
which:

- injects `0.11.0` (the tag, minus any leading `v`) into `Cargo.toml`,
- builds wheels for Linux (glibc + musl, x86_64/aarch64), macOS
  (x86_64/arm64), and Windows (x64/arm64), plus an sdist,
- runs the test suite against each wheel,
- publishes everything to PyPI.

Platforms without a prebuilt wheel (e.g. 32-bit, s390x/ppc64le) install from the
sdist, which compiles the Rust extension and therefore needs a Rust toolchain.

`pip install seqfold` picks up the new version a few minutes after the run
finishes.

### Publishing auth: PyPI Trusted Publishing (no API token)

Publishing uses [PyPI Trusted Publishing](https://docs.pypi.org/trusted-publishers/)
over GitHub OIDC — there is **no API token to store or rotate**. Two one-time
setup steps:

1. **GitHub environment.** Create an Environment named **`pypi`** under
   repo *Settings → Environments*. The `release` job runs in it
   (`environment: pypi`). Add **required reviewers** there if you want a manual
   approval gate before anything reaches PyPI — useful when maintainers have
   commit access but shouldn't all have publish rights.

2. **PyPI trusted publisher.** On PyPI: *Project → Settings → Publishing → Add a
   GitHub publisher*, with:

   | Field | Value |
   | --- | --- |
   | Owner | `Lattice-Automation` |
   | Repository | `seqfold` |
   | Workflow name | `release.yml` |
   | Environment name | `pypi` |

After that, any published GitHub Release builds and publishes automatically. (If
you previously added a `PYPI_API_TOKEN` repo secret for the old token-based flow,
you can delete it — it's no longer used.)

### Dry run

Trigger the workflow manually (*Actions → Release → Run workflow*, i.e.
`workflow_dispatch`) to build the full wheel matrix without publishing — the
publish step only runs on a `release` event.

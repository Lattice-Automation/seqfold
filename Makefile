.PHONY: examples test install develop build profile parse

PY ?= python3

# Build the Rust extension and install it (editable) into the active env.
develop:
	maturin develop --release

# Build a release wheel under target/wheels.
build:
	maturin build --release

# Build + install from source (pip uses the maturin backend in pyproject.toml).
install:
	$(PY) -m pip install .

# Rust unit tests (ported internal tests) + Python public-API tests.
test:
	cargo test --lib
	$(PY) -m unittest discover tests -p '*_test.py'

# Regenerate src/core/data.rs from the original Python energy tables.
# Requires the original seqfold/dna.py and seqfold/rna.py (see git history).
parse:
	$(PY) ./codegen/gen_data.py

.PHONY: examples
examples: develop
	$(PY) ./tests/fold_examples.py

profile: develop
	$(PY) ./tests/fold_profile.py

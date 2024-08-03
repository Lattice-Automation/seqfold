.PHONY: examples

PY ?= python3

install:
	rm -f ./seqfold/*.c
	$(PY) -m pip install .

test:
	$(PY) -m unittest discover tests -p '*_test.py'

parse:
	$(PY) ./data/_rna.py
	black ./seqfold/rna.py

patch: test
	bumpversion patch
	$(PY) setup.py sdist bdist_wheel
	$(PY) -m twine upload dist/* --skip-existing

minor: test
	bumpversion minor
	$(PY) setup.py sdist bdist_wheel
	$(PY) -m twine upload dist/* --skip-existing

.PHONY: examples
examples: install
	$(PY) ./tests/fold_examples.py

profile: install
	$(PY) ./tests/fold_profile.py

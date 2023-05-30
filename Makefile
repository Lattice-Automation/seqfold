.PHONY: examples

PY ?= python3

install:
	rm -f ./seqfold/*.c
	$(PY) setup.py install

test:
	$(PY) -m unittest discover tests -p '*_test.py'

parse:
	$(PY) ./data/_rna.py
	black ./seqfold/rna.py

patch: test
	bumpversion patch
	python3 setup.py sdist bdist_wheel
	python3 -m twine upload dist/* --skip-existing

minor: test
	bumpversion minor
	python3 setup.py sdist bdist_wheel
	python3 -m twine upload dist/* --skip-existing

examples: install
	$(PY) ./tests/fold_examples.py

profile: install
	$(PY) ./tests/fold_profile.py

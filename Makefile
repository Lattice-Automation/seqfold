.PHONY: examples

install:
	python3 setup.py install

test:
	python3 -m unittest discover tests -p '*_test.py'

parse:
	python3 ./data/_rna.py
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
	python3 ./tests/fold_examples.py

profile: install
	python3 ./tests/fold_profile.py

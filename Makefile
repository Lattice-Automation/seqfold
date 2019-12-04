install:
	python setup.py install
	pip3 install . --user

test:
	python3 -m unittest discover tests -p '*_test.py'

parse:
	python3 ./scripts/rna.py
	black ./seqfold/rna.py

patch: test
	bumpversion patch
	python3 setup.py sdist bdist_wheel
	python3 -m twine upload dist/* --skip-existing

minor: test
	bumpversion minor
	python3 setup.py sdist bdist_wheel
	python3 -m twine upload dist/* --skip-existing


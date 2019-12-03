install:
	python setup.py install
	pip3 install . --user

test:
	python3 -m unittest discover tests -p '*_test.py'

patch: test
	bumpversion patch

minor: test
	bumpversion minor


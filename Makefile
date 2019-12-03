install:
	python setup.py install
	pip3 install . --user

test:
	python setup.py test
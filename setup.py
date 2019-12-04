import sys

from pkg_resources import VersionConflict, require
from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("requirements.txt") as f:
    requirements = f.read().splitlines()

try:
    require("setuptools>=38.3")
except VersionConflict:
    print("Error: version of setuptools is too old (<38.3)!")
    sys.exit(1)

setup(
    name="seqfold",
    version="0.3.0",
    description="Predict the minimum free energy structure of nucleic acids",
    author="JJTimmons",
    author_email="jtimmons@latticeautomation.com",
    license="mit",
    packages=find_packages(exclude=["data"]),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Lattice-Automation/seqfold",
    install_requires=requirements,
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    entry_points={"console_scripts": ["seqfold=seqfold.fold:run"],},
)

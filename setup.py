#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""
from setuptools import setup, find_packages

import versioneer

with open("docs/readme.rst") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()

with open("requirements.txt") as requirements_file:
    requirements = requirements_file.read().split("\n")
    if "" in requirements:
        requirements.remove("")


setup_requirements = [
    "pytest-runner",
]

test_requirements = ["pytest"]
extras_require = {"plots": ["rdkit-pypi==2021.3.4"]}

setup(
    author="Data Revenue GmbH",
    author_email="carlos@datarevenue.com",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ],
    description="A package to access Omigami services.",
    entry_points={
        "console_scripts": [
            "omigami=omigami.cli:omigami",
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + "\n\n" + history,
    include_package_data=True,
    keywords="omigami",
    name="omigami-client",
    packages=find_packages(include=["omigami"]),
    setup_requires=setup_requirements,
    test_suite="tests",
    tests_require=test_requirements,
    url="https://github.com/omigami/omigami",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    zip_safe=False,
    extras_require=extras_require,
)

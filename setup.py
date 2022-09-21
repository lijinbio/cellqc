#!/usr/bin/env python3

import os
import io
import re
import ast

from setuptools import setup, find_packages

DEPENDENCIES = ['click', 'snakemake', 'pygraphviz']
EXCLUDE_FROM_PACKAGES = ["contrib", "docs", "tests*"]
CURDIR = os.path.abspath(os.path.dirname(__file__))

with io.open(os.path.join(CURDIR, "README.md"), "r", encoding="utf-8") as readme_file:
    readme = readme_file.read()

def get_version():
    main_file = os.path.join(CURDIR, "cellqc", "cellqc.py")
    _version_re = re.compile(r"__version__\s+=\s+(?P<version>.*)")
    with open(main_file, "r", encoding="utf8") as f:
        match = _version_re.search(f.read())
        version = match.group("version") if match is not None else '"unknown"'
    return str(ast.literal_eval(version))

setup(
    name="cellqc",
    author="Jin Li",
    author_email="lijin.abc@gmail.com",
    python_requires=">=3.6",
    description="Cellqc standardizes the qualiy control of single-cell RNA-Seq (scRNA) data to render clean feature count matrices.",
    install_requires=DEPENDENCIES,
    packages=find_packages(exclude=EXCLUDE_FROM_PACKAGES),
    long_description=readme,
    long_description_content_type="text/markdown",
    include_package_data=True,
    keywords=[],
    scripts=[],
    setup_requires=[],
    entry_points={"console_scripts": ["cellqc=cellqc.cellqc:main"]},
    url="https://github.com/lijinbio/cellqc",
    version=get_version(),
    zip_safe=False,
    license="MIT license",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
    ],
)

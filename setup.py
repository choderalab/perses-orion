# (C) 2018-2019 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.

from re import compile
from ast import literal_eval
from sys import argv, exit
from json import dumps
from setuptools import setup, find_packages, convert_path

from pip._internal.req import parse_requirements
from pip._internal.download import PipSession
def get_reqs(reqs):
  return [str(ir.req) for ir in reqs]
requirements = get_reqs(parse_requirements("orion-requirements.txt", session=PipSession()))

if argv[-1] == "--requires":
    print(dumps(requirements))
    exit()

# Obtain version of the package
_version_re = compile(r"__version__\s+=\s+(.*)")
version_file = convert_path("./cubes/__init__.py")
with open(version_file, "rb") as f:
    version = str(literal_eval(_version_re.search(f.read().decode("utf-8")).group(1)))


setup(
    name="perses-orion",
    version=version,
    packages=find_packages(exclude=["tests/*", "floes/*"]),
    author="John D. Chodera",
    author_email="john.chodera@choderalab.org",
    description="Relative alchemical free energy calculations with perses on Orion",
    license="Other/Proprietary License",
    keywords="openeye cloud orion",
    include_package_data=True,
    install_requires=requirements,
    classifiers=[
        "Development Status :: Beta",
        "Intended Audience :: Orion",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python",
        "Topic :: Software Development :: Libraries",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.6",
    ],
)

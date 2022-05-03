#!/usr/bin/python
# -*- coding: utf-8 -*-

from setuptools import find_packages, setup

def parse_requirements(filename):
    lineiter = (line.strip() for line in open(filename))
    return [line for line in lineiter if line and not line.startswith("#")]

setup(
    name="Cals",
    version="1.0.0",
    description='Calculate stellar chromospheric activity parameters and generate corresponding spectrum diagram of Ca II H&K lines with LAMOST LRS spectra',
    license = 'BSD-3-Clause License',
    #install_requires=parse_requirements('requirements.txt'),
    packages=find_packages()
)
    
    

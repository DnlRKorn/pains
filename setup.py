#!/usr/bin/env python

import sys
from setuptools import setup

setup(name='pains',
      version='0.1',
      install_requires=[
          "flask",
		  "flask_restful",
		  "rdkit"		  
      ],
      packages=["pains"]
      )

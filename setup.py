#!/usr/bin/env python3

from setuptools import setup, find_packages


setup(name='ars_rqc',
      version='0.1.0.dev1',
      description='A rolling quality control workflow for the USDA Agricultural \
      Research Service, Genomics and Bioinformatics Research Unit',
      url='https://github.com/arivers/ars_rqc',
      author='Adam R. Rivers',
      author_email='adam.rivers@ars.usda.gov',
      license='CC0 1.0',
      packages=find_packages(),
      install_requires=[
          'pandas',
          'numpy'
          ],
      test_suite='nose.collector',
      tests_require=['nose'],
      scripts=['bin/rqcfilter.py'],
      zip_safe=False)

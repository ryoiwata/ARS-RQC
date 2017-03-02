#!/usr/bin/env python3

from setuptools import setup


setup(name='ars_rqc',
      version='0.1',
      description='A rolling quality control worflow for the USDA Agricultural \
      Research Service, Geomics and Bioinformatics Research Unit',
      url='https://github.com/arivers/ars_rqc',
      author='Adam R. Rivers',
      author_email='adam.rivers@ars.usda.gov',
      license='CC0 1.0',
      packages=['ars_rqc'],
      install_requires=[
          'pandas'
          ],
      test_suite='nose.collector',
      tests_require=['nose'],
      scripts=['bin/rqcfilter.py'],
      zip_safe=False)

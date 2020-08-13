#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst


import sys
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read().replace('.. :changelog:', '')

with open('requirements.txt') as requirements_file:
    requirements = requirements_file.read()


setup(
    name='tripexPro',
    version='0.1',
    description='It processes the data vertically pointed radars',
    long_description=readme + '\n\n' + history,
    author='Jose Dias Neto',
    author_email='jdiasn@gmail.com',
    url='http://jdiasn.github.io',
    packages=[
        'tripexPro',
    ],
    package_dir = {'tripexPro':
                   'tripexPro'},
    license='3-clause BSD',
    install_requires=requirements,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'License :: OSI Approved :: BSD License',
        ],
    keywords='',
    include_package_data=True,
    zip_safe=False,
    extras_require = {
        'plot': ["matplotlib"]
    }
)

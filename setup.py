#!/usr/bin/env python

from setuptools import setup

setup(name='hicplotlib',
      version='0.2.2',
      description='Hi-C data plotting and analysis tool',
      author='Ilya Flyamer',
      author_email='flyamer@gmail.com',
      url='https://github.com/Phlya/hicplotlib',
      packages=['hicplotlib'],
      install_requires=['numpy', 'pandas', 'matplotlib', 'seaborn', 'pybedtools', 'scipy']
     )

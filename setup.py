"""
Nothing too special happening here, typical setup.py

Main thing is to define the entrypoints for the GUI and CLI.
"""

from setuptools import setup, find_packages

with open("asarix/__init__.py") as f:
    exec([x for x in f.readlines() if '__version__' in x][0])

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("requirements.txt", "r") as f:
    requirements = f.read()

setup(
  name='asarix',
  version=__version__,
  author='Shuzhao Li',
  author_email='shuzhao.li@gmail.com',

  description='eXposome data mining based on asari',
  long_description=long_description,
  long_description_content_type="text/markdown",
  url='https://github.com/shuzhao-li/asari-x',
  license='BSD 3-Clause',
  keywords='exposome metabolomics bioinformatics mass spectrometry',

  classifiers=[
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Natural Language :: English',
    'Operating System :: OS Independent',
    'Programming Language :: Python :: 3',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Scientific/Engineering :: Chemistry',
    'Topic :: Software Development :: Libraries :: Python Modules',
  ],

  packages=find_packages(
    include=['*', '']
  ),
  # data_files=[ ('asari/db', ['asari/db/mass_indexed_compounds.pickle', 'asari/db/emp_cpds_trees.pickle']) ],
  include_package_data=True,
  zip_safe=True,
  entry_points = {
        'console_scripts': ['asarix=asarix.main:cli', 'asarix-gui=asarix.gui:main_gui'],
    },

  python_requires='>=3.7',
  install_requires=requirements,

)

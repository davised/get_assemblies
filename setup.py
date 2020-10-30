from setuptools import find_packages, setup

with open('get_assemblies/__version__.py', 'r') as f:
    for line in f:
        if line.startswith('__version__'):
            version = line.strip().split('=')[1].strip(' \'"')
            break
    else:
        version = '0.0.1'

with open('README.rst', 'r', encoding='utf-8') as f:
    readme = f.read()

REQUIRES = ['rich']

kwargs = {
    'name': 'get_assemblies',
    'version': version,
    'description': 'NCBI E-utilities wrapper for assembly downloads',
    'long_description_content_type': 'text/x-rst',
    'long_description': readme,
    'author': 'Edward Davis',
    'author_email': 'ed@cgrb.oregonstate.edu',
    'maintainer': 'Edward Davis',
    'maintainer_email': 'ed@cgrb.oregonstate.edu',
    'url': 'https://github.com/davised/get_assemblies',
    'license': 'Custom',
    'classifiers': [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'License :: Free for non-commercial use',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    'install_requires': REQUIRES,
    'tests_require': ['coverage', 'pytest'],
    'packages': find_packages(exclude=('tests', 'tests.*')),

}

setup(**kwargs,
      entry_points={
          "console_scripts": ['get_assemblies = get_assemblies.__main__:main']
      }
      )

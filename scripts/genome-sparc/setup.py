from setuptools import setup, find_packages
import os.path

# Get the version:
version = {}
with open(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'safeports', 'version.py')) as f: exec(f.read(), version)

setup(
    name = 'genome-sparc',
    version = version['__version__'],
    description = 'Annotate genomic regions that are safe for gene editing and recode the region.',
    author = 'Juliane Weller',
    author_email = 'jw38@sanger.ac.uk',
    url = 'https://github.com/UPDATE_THIS',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 3'
    ],
    packages = find_packages(),
    install_requires = [
    ],
    python_requires = '>=3',
    entry_points = {
        'console_scripts': [
            'sparc=genome-sparc.scripts:main'
        ]
    }
)

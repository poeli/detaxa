from setuptools import setup, find_packages

from src.pytaxa import __version__ as ptv

INSTALL_REQUIRES = [
    pkg for pkg in open('requirements.txt').readlines()
]

setup(
    name='pytaxa',
    version= ptv,
    author='Paul Li',
    author_email='po-e@lanl.gov',
    url='https://github.com/poeli/taxonomy',
    license='LICENSE',
    packages=['pytaxa'],
    # packages=find_packages(include=['pytaxa', 'pytaxa.*']),
    include_package_data=True,
    package_dir={'': 'src'},
    package_data={'pytaxa':['taxonomy_db/*.json']},
    description='NCBI taxonomy and lineage lookup',
    keywords=["taxonomy"],
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    entry_points={ 
        'console_scripts': ['pytaxa = pytaxa.__main__:cli' ] 
    },
    install_requires=INSTALL_REQUIRES,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Operating System :: OS Independent",
        "Development Status :: 4 - Beta",
    ],
    zip_safe=False
)
#!/usr/bin/env python3
"""
FoTChemistry Setup Script

Setup script for the Field of Truth Chemistry project.
Ensures deterministic installation and dependency management.
"""

from setuptools import setup, find_packages
import os
import sys

# Ensure Python 3.9+
if sys.version_info < (3, 9):
    raise RuntimeError("FoTChemistry requires Python 3.9 or higher")

# Read long description from README
def read_file(filename):
    here = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(here, filename), 'r', encoding='utf-8') as f:
        return f.read()

# Read requirements
def read_requirements(filename):
    requirements = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                requirements.append(line)
    return requirements

# Package metadata
NAME = "fot-chemistry"
VERSION = "0.1.0"
DESCRIPTION = "Field of Truth open lab notebook and truth ledger for chemistry"
LONG_DESCRIPTION = read_file("README.md") if os.path.exists("README.md") else DESCRIPTION
AUTHOR = "FoT Research Team"
AUTHOR_EMAIL = "research@fieldoftruth.org"
URL = "https://github.com/FortressAI/FoTChemistry"
LICENSE = "Apache-2.0"

# Package requirements
INSTALL_REQUIRES = read_requirements("requirements.txt")

# Development requirements
EXTRAS_REQUIRE = {
    'dev': [
        'pytest>=7.4.3',
        'pytest-cov>=4.1.0',
        'pytest-benchmark>=4.0.0',
        'black>=23.11.0',
        'flake8>=6.1.0',
        'mypy>=1.7.1',
        'pre-commit>=3.5.0',
    ],
    'docs': [
        'sphinx>=7.2.6',
        'sphinx-rtd-theme>=1.3.0',
        'myst-parser>=2.0.0',
        'mkdocs>=1.5.3',
        'mkdocs-material>=9.4.8',
    ],
    'jupyter': [
        'jupyter>=1.0.0',
        'jupyterlab>=4.0.8',
        'ipywidgets>=8.1.1',
    ],
    'visualization': [
        'matplotlib>=3.8.1',
        'seaborn>=0.13.0',
        'plotly>=5.17.0',
        'nglview>=3.0.8',
        'py3dmol>=2.0.4',
    ],
    'chemistry': [
        'rdkit>=2023.9.1',
        'psi4>=1.7.0',
        'openmm>=8.0.0',
        'qcengine>=0.27.0',
        'cclib>=1.8.0',
        'openbabel>=3.1.1',
    ],
    'ml': [
        'torch>=2.1.0',
        'pytorch-lightning>=2.1.0',
        'dgl>=1.1.2',
        'torch-geometric>=2.4.0',
        'moleculenet>=0.1.0',
    ],
    'graph': [
        'neo4j>=5.13.0',
        'rdflib>=7.0.0',
        'sparqlwrapper>=2.0.0',
        'owlready2>=0.45',
    ],
    'pipelines': [
        'snakemake>=7.32.4',
        'cwltool>=3.1.20231013155012',
        'docker>=6.1.3',
        'apptainer>=1.2.4',
    ]
}

# All optional dependencies
EXTRAS_REQUIRE['all'] = [
    item for sublist in EXTRAS_REQUIRE.values() 
    for item in sublist
]

# Package classifiers
CLASSIFIERS = [
    "Development Status :: 3 - Alpha",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: Apache Software License",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "Topic :: Scientific/Engineering :: Chemistry",
    "Topic :: Scientific/Engineering :: Artificial Intelligence",
    "Topic :: Software Development :: Libraries :: Python Modules",
    "Topic :: Database :: Database Engines/Servers",
    "Operating System :: MacOS :: MacOS X",
    "Operating System :: POSIX :: Linux",
    "Operating System :: Microsoft :: Windows",
]

# Keywords
KEYWORDS = [
    "chemistry",
    "knowledge graph", 
    "field of truth",
    "computational chemistry",
    "cheminformatics",
    "agentic systems",
    "truth mining",
    "ontology",
    "FAIR data",
    "reproducible research",
    "open science",
    "neo4j",
    "RDKit",
    "quantum chemistry"
]

# Python requirements
PYTHON_REQUIRES = ">=3.9"

# Package data
PACKAGE_DATA = {
    'fot': [
        'ontology/*.owl',
        'ontology/*.ttl',
        'ontology/*.jsonld',
        'data/*.json',
        'data/*.csv',
        'configs/*.yaml',
        'sparql/*.rq',
    ]
}

# Data files
DATA_FILES = [
    ('ontology', ['ontology/fot_chemistry.owl']),
]

# Entry points for command line tools
ENTRY_POINTS = {
    'console_scripts': [
        'fot-chem=fot.cli:main',
        'fot-setup-graph=fot.setup:setup_graph_db',
        'fot-validate-claims=fot.validate:main',
        'fot-mine-truth=fot.truth_mining:main',
        'fot-export-zenodo=fot.export:zenodo_export',
    ],
}

# Custom commands
class DeterministicInstallCommand:
    """Custom installation command that ensures deterministic setup"""
    
    def run(self):
        import subprocess
        import sys
        
        # Verify RDKit availability
        try:
            from rdkit import Chem
            print("✓ RDKit chemistry toolkit available")
        except ImportError:
            print("⚠ Warning: RDKit not yet installed - essential for chemistry operations")
        
        # Verify Neo4j connectivity (optional)
        try:
            from neo4j import GraphDatabase
            print("✓ Neo4j graph database driver available")
        except ImportError:
            print("⚠ Warning: Neo4j driver not available - required for knowledge graph")
        
        # Set deterministic environment variables
        os.environ['PYTHONHASHSEED'] = '1337'
        
        print("✓ Deterministic environment configured")
        print("✓ FoTChemistry installation completed")

def setup_package():
    """Main setup function"""
    
    # Verify system requirements
    print("Setting up FoTChemistry...")
    print(f"Python version: {sys.version}")
    print(f"Platform: {sys.platform}")
    
    # Mac M4 specific checks
    if sys.platform == "darwin":
        import platform
        machine = platform.machine()
        if machine == "arm64":
            print("✓ Apple Silicon (M-series) detected - optimized chemistry libraries available")
        else:
            print(f"⚠ Warning: Expected arm64, got {machine}")
    
    setup(
        name=NAME,
        version=VERSION,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        long_description_content_type="text/markdown",
        author=AUTHOR,
        author_email=AUTHOR_EMAIL,
        url=URL,
        license=LICENSE,
        packages=find_packages(exclude=['tests*', 'docs*', 'examples*']),
        package_data=PACKAGE_DATA,
        data_files=DATA_FILES,
        include_package_data=True,
        install_requires=INSTALL_REQUIRES,
        extras_require=EXTRAS_REQUIRE,
        python_requires=PYTHON_REQUIRES,
        classifiers=CLASSIFIERS,
        keywords=' '.join(KEYWORDS),
        entry_points=ENTRY_POINTS,
        zip_safe=False,  # Required for package data access
        
        # Custom metadata
        project_urls={
            "Bug Reports": f"{URL}/issues",
            "Source": URL,
            "Documentation": f"{URL}/docs",
            "Funding": "https://github.com/sponsors/fot-research",
            "Wiki": f"{URL}/wiki",
        },
        
        # Setuptools options
        options={
            'build_py': {
                'compile': True,
                'optimize': 2,
            },
        },
        
        # Additional metadata for PyPI
        maintainer=AUTHOR,
        maintainer_email=AUTHOR_EMAIL,
        platforms=["macOS", "Linux", "Windows"],
    )

if __name__ == "__main__":
    setup_package()
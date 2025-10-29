# setup.py
from setuptools import setup, find_packages

# If you keep a README.md, this will use it as the long description.
def read_long_description():
    try:
        with open("README.md", encoding="utf-8") as f:
            return f.read()
    except FileNotFoundError:
        return ""

setup(
    name="GeneScanner",        
    version="1.0.0",                 
    description="Mutation analysis on nucleotide/protein alignments with Excel reports",
    long_description=read_long_description(),
    long_description_content_type="text/markdown",
    author="Seungwon Ko",
    author_email="kell7366@ox.ac.uk",
    url="https://github.com/jeju2486/GeneScanner",
    license="MIT",                 
    license_files=("LICENSE",),              
    python_requires=">=3.10",                 
    packages=find_packages(include=('GeneScanner', 'GeneScanner.*')),
    install_requires=[
        "biopython>=1.81",
        "numpy>=1.23",
        "pandas>=1.5",
        "XlsxWriter>=3.0",
    ],
    entry_points={
        "console_scripts": [
            "GeneScanner = GeneScanner.mutation_analysis:main",
            "GeneScanner-viewer = GeneScanner.genescanner_viewer:main"
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Environment :: Console",
    ],
    project_urls={
        "Source": "https://github.com/jeju2486/GeneScanner",
        "Tracker": "https://github.com/jeju2486/GeneScanner/issues",
    },
    include_package_data=True,
    zip_safe=False,
)


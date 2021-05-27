from setuptools import setup, Extension, find_packages
import numpy
from Cython.Build import cythonize

import glob

setup(
    name='reat',
    version='0.2.0',
    packages=find_packages('.', exclude=["tests"]),
    url='https://github.com/ei-corebioinformatics/reat',
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific Engineering :: Bio/Informatics",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8"
    ],
    license='MIT',
    author='Luis Yanes',
    author_email='luis.yanes@earlham.ac.uk',
    description='Robust Eukaryotic Annotation Toolkit',
    zip_safe=False,
    keywords="gene annotation WDL pipeline workflow cromwell transcriptome homology",
    ext_modules=cythonize(
        [
            Extension("annotation.lib.cy_utils.contrast",
                      include_dirs=[numpy.get_include()],
                      sources=["annotation/lib/cy_utils/contrast.pyx"]),
        ],
        compiler_directives={"language_level": "3"}
    ),
    scripts=[
        script for script in glob.glob("annotation/scripts/*")
    ],
    install_requires=[
        "biopython~=1.78",
        "mikado~=2.3.0",
        "pyfaidx~=0.5.8",
        "numpy~=1.20.3",
        "jsonschema~=3.2.0",
        "pyyaml~=5.4.1",
        "parasail~=1.2.4",
        '2passtools @ git+https://github.com/bartongroup/2passtools.git#d4378d0'
    ],
    extras_require={
        'docs': [
            'sphinx',
            'sphinx_pdj_theme',
        ]
    },
    package_data={
        "validation": ["transcriptome.schema.json", "homology.schema.json"],
        "annotation": ["transcriptome_module/*.wdl",
                       "transcriptome_module/*/**/*.wdl",
                       "transcriptome_module/*/**/**/*.wdl",
                       "homology_module/*.wdl"],
        "annotation.lib.cy_utils": ["annotation/lib/cy_utils/contrast.pxd"]
    },
    entry_points={
        "console_scripts": [
            "reat=annotation.__main__:main"
        ]
    }
)

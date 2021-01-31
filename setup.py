from setuptools import setup, find_packages
import glob

setup(
    name='reat',
    version='0.0.7',
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
    keywords="gene annotation WDL pipeline workflow",
    scripts=[
        script for script in glob.glob("annotation/scripts/*")
    ],
    install_requires=[
        "pyyaml~=5.3",
        "pyfaidx~=0.5.8",
        "jsonschema~=3.2.0",
        "biopython~=1.78"
    ],
    package_data={
        "validation": ["transcriptome.schema", "homology.schema"],
        "annotation": ["transcriptome_module/*.wdl",
                       "transcriptome_module/*/**/*.wdl",
                       "transcriptome_module/*/**/**/*.wdl",
                       "homology_module/*.wdl"],
    },
    entry_points={
        "console_scripts": [
            "reat=annotation.__main__:main"
        ]
    }
)

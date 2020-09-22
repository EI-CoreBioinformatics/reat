from setuptools import setup, find_packages
import glob

setup(
    name='ei-cautious-broccoli',
    version='0.2',
    packages=find_packages('./'),
    url='github.com/ei-corebioinformatics/ei-cautious-broccoli',
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific Engineering :: Bio/Informatics",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        'Programming Language :: Python :: 3.4',
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8"
    ],
    license='MIT',
    author='Luis Yanes',
    author_email='luis.yanes@earlham.ac.uk',
    description='Robust and extendable eukaryotic annotation toolkit',
    zip_safe=False,
    keywords="gene annotation WDL pipeline workflow",
    scripts=[
        script for script in glob.glob("annotation/scripts/*")
    ],
    install_requires=[
        "pyyaml~=5.3",
        "pyfaidx~=0.5.8",
        "jsonschema~=3.2.0",
    ],
    package_data={
        "validation": ["reat.schema"],
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

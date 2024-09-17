from setuptools import setup, find_packages, Extension


NAME = "DEPlots"
DESCRIPTION = "Package for identification of RNA dependent Proteins from mass spec data "

LONGDESC = DESCRIPTION #Todo: write Readme


setup(
    name=NAME,
    version="0.1.1",
    author="domonik",
    author_email="dominik.rabsch@gmail.com",
    packages=find_packages(),
    package_dir={"DEplots": "./DEplots"},
    license="LICENSE",
    url="https://github.com/domonik/DEPlots",
    description=DESCRIPTION,
    long_description=LONGDESC,
    long_description_content_type="text/markdown",
    include_package_data=True,
    package_data={
        "DEplots.dashboard": ["assets/*", "default_config.yaml"],
    },
    install_requires=[
        "numpy",
        "plotly>=5.16",
        "pandas",
        "dash>=2.5",
        "dash_bootstrap_components",
        "pyYAML",
        "gffutils",
        "pysam"
    ],
    scripts=[
        "DEplots/executables.py",
    ],
    entry_points={
        "console_scripts": [
            "DEPlots = DEplots.executables:main"
        ]
    },
)
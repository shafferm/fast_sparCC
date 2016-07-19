from setuptools import setup, find_packages

__author__ = 'shafferm'


setup(
    name="sparcc_fast",
    version="0.1",
    install_requires=["numpy", "scipy", "biom-format", "pandas"],
    scripts=["scripts/sparCC_fast.py"],
    packages=find_packages()
)
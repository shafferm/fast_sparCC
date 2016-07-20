from setuptools import setup, find_packages

__author__ = 'shafferm'


setup(
    name="sparcc_fast",
    version="0.1.1",
    author="Michael Shaffer",
    author_email="michael.shaffer@ucdenver.edu",
    description="A fast command line interface to find correlations in biom tables with SparCC.",
    license="BSD",
    url="https://github.com/shafferm/fast_sparCC",

    install_requires=["numpy", "scipy", "biom-format", "pandas"],
    scripts=["scripts/fast_sparCC.py"],
    packages=find_packages()
)
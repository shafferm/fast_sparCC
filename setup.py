from setuptools import setup, find_packages

__author__ = 'shafferm'


setup(
    name="fast_sparCC",
    version="v0.1.6",
    author="Michael Shaffer",
    author_email="michael.shaffer@ucdenver.edu",
    description="A fast command line interface to find correlations in biom tables with SparCC.",
    license="BSD",
    url="https://github.com/shafferm/fast_sparCC",
    download_url="https://github.com/shafferm/fast_sparCC/archive/v0.1.6.tar.gz",

    install_requires=["numpy", "scipy", "biom-format", "pandas"],
    scripts=["scripts/fast_sparCC.py"],
    packages=find_packages()
)

from setuptools import setup, find_packages

VERSION = '0.1.0.dev0'

with open('README.rst', 'r') as readme:
    README_TEXT = readme.read()

with open("requirements.txt") as f:
    INSTALL_REQUIRES = f.readlines()

# main setup configuration class
setup(
    name='SimPhoNy Tools Collection',
    version=VERSION,
    author='SimPhoNy, EU FP7 Project (Nr. 604005) www.simphony-project.eu',
    description='Collection of tools for SimPhoNy',
    long_description=README_TEXT,
    install_requires=INSTALL_REQUIRES,
    packages=find_packages(),
    entry_points={'simphony': ['tools = tools']}, 
)

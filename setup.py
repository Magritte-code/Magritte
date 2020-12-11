import os
import re
import sys
import shutil

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext


with open("README.md", "r") as file:
    long_description = file.read()


with open("src/configure.hpp", "r") as file:
    for line in file.readlines():
        if 'MAGRITTE_VERSION' in line:
            # Get the version, which is between ""
            __version__ = re.findall('"(.*?)"', line)[0]


class my_build_ext(build_ext):
    def build_extension(self, ext):
        '''
        Copy the already-compiled core.so file!
        '''
        shutil.copyfile('bin/core.so', self.get_ext_fullpath(ext.name))


setup(
    name         ="magritte",
    version      =__version__,
    author       ="Frederik De Ceuster",
    author_email ="frederik.deceuster@gmail.com",
    description="A modern software library for simulating radiation transport.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    ext_modules=[Extension('magritte/core', sources=[])],
    packages=find_packages(),
    install_requires=[
        'h5py>=2.9',
        'numpy>=1.19',
        'scipy>=1.3',
        'astropy>=4.0',
        'numba>=0.46',
        'meshio>=4.3',
        'healpy>=1.13',
        'jupyterlab>=2.2',
        'yt>=3.4',
        'mpi4py>=3',
        'tqdm>=4.42',
    ],
    cmdclass={'build_ext': my_build_ext},
    url="https://github.com/Magritte-code/Magritte",
    zip_safe=False,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Environment :: GPU",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: C++",
        "Programming Language :: Python :: 3",
        "Natural Language :: English",
        "Operating System :: Unix",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics"
    ],
    python_requires='>=3.6',
)

import os
import re
import sys
import shutil

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from subprocess import call

this_dir = os.path.dirname(os.path.realpath(__file__))
print(os.listdir(os.getcwd()))

with open(os.path.join(this_dir, "./README.md"), "r") as file:
    long_description = file.read()

with open("./src/configure.hpp", "r") as file:
    for line in file.readlines():
        if 'MAGRITTE_VERSION' in line:
            # Get the version, which is between quotes ("")
            __version__ = re.findall('"(.*?)"', line)[0]


class my_build_ext(build_ext):

    def build_extension(self, ext):
        '''
        Build magritte C++ core module.
        '''
        # Compile module
        call('bash ./compile.sh', shell=True)
        # Copy the already-compiled core.so file.
        shutil.copyfile('./bin/core.so'), self.get_ext_fullpath(ext.name))


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
    include_package_data=True,
    install_requires=[
        'h5py',
        'numpy',
        'scipy',
        'astropy',
        'numba',
        'meshio',
        'healpy',
        'yt',
        'mpi4py',
        'tqdm',
        'ipywidgets',
        'palettable',
        'plotly'
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
    python_requires='>=3.7',
)

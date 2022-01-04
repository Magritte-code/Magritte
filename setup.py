import os
import re
import sys
import shutil

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from subprocess import call

this_dir = os.path.dirname(os.path.realpath(__file__))

with open(os.path.join(this_dir, "README.md"), "r") as file:
    long_description = file.read()


with open(os.path.join(this_dir, "CMakeLists.txt"), "r") as file:
    # Extract the lines containing the project description
    lines = re.findall('project.*\([^\)]*\)', file.read())[0]
    # Extract the version nuber form those lines
    __version__ = re.findall('\d*\.\d*\.\d*', lines)[0]


class my_build_ext(build_ext):

    def build_extension(self, ext):
        '''
        Build magritte C++ core module.
        '''
        # Compile module
        call(f'bash {os.path.join(this_dir, "compile.sh")}', shell=True)
        # Copy the already-compiled core.so file.
        shutil.copyfile(os.path.join(this_dir, 'bin', 'core.so'), self.get_ext_fullpath(ext.name))


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

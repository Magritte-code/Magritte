# Used example:
# https://github.com/sizmailov/pyxmolpp2/blob/master/setup.py

with open("README.md", "r") as file:
    long_description = file.read()

# with open("CMakeLists.txt", "r") as file:
    # long_description = file.read()


__version__ = '0.0.1'


import os
import re
import sys
import platform
import subprocess

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion

#class CoreExtension(Extension):
#    def __init__(self, name, sourcedir=''):
#        Extension.__init__(self, name, sources=[])
#        self.sourcedir = os.path.abspath(sourcedir)
#
# class CMakeExtension(Extension):
#     def __init__(self, name, sourcedir=''):
#         Extension.__init__(self, name, sources=[])
#         self.sourcedir = os.path.abspath(sourcedir)
#
#
# class CMakeBuild(build_ext):
#
#     def run(self):
#         # Check if CMake is installed
#         try:
#             out = subprocess.check_output(['cmake', '--version'])
#         except OSError:
#             raise RuntimeError(
#                 "CMake must be installed to build the following extensions: "
#                 + ", ".join(e.name for e in self.extensions)
#             )
#         # Check CMake version
#         cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
#         if cmake_version < LooseVersion('3.18.0'):
#             raise RuntimeError("CMake >= 3.18.0 is required")
#         # Build extension
#         for ext in self.extensions:
#             self.build_extension(ext)
#
#     def build_extension(self, ext):
#         extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
#         cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
#                       '-DOMP_PARALLEL=ON',
#                       '-DMPI_PARALLEL=OFF',
#                       '-DGPU_ACCELERATION=OFF',
#                       '-DPYTHON_EXECUTABLE:FILEPATH=' + sys.executable]
#
#         # build_type = os.environ.get("BUILD_TYPE", "Release")
#         # build_args = ['--config', build_type]
#
#         # Pile all .so in one place and use $ORIGIN as RPATH
#         # cmake_args += ["-DCMAKE_BUILD_WITH_INSTALL_RPATH=TRUE"]
#         # cmake_args += ["-DCMAKE_INSTALL_RPATH={}".format("$ORIGIN")]
#
#         # cmake_args += ['-DCMAKE_BUILD_TYPE=' + build_type]
#         # build_args += ['--', '-j4']
#
#         # env = os.environ.copy()
#         # env['CXXFLAGS'] = '{}'.format(env.get('CXXFLAGS', ''))
#         # if not os.path.exists(self.build_temp):
#             # os.makedirs(self.build_temp)
#         # subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
#         # subprocess.check_call(['cmake',
#                                # '--build', '../',
#                                # '--target', 'core'#ext.name
#                                # ] + build_args,
#                               # cwd=self.build_temp)

class my_build_ext(build_ext):
    def build_extension(self, ext):
        '''
        Copy the already-compiled core.so file!
        '''
        import shutil
        import os.path
        print('HERE')
        # os.makedirs(os.path.dirname(self.get_ext_fullpath(ext.name)))
        print('HERE2')
        print(ext.name)

        shutil.copyfile('bin/core.so', self.get_ext_fullpath(ext.name))
        print('HERE3')



setup(
    name         ="magritte",
    version      =__version__,
    author       ="Frederik De Ceuster",
    author_email ="frederik.deceuster@gmail.com",
    description="A modern software library for simulating radiation transport.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    # ext_modules=[CMakeExtension('magritte')],
    ext_modules=[Extension('magritte/core', sources=[])],
    packages=find_packages(),
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

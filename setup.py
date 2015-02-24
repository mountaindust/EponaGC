# -*- coding: utf-8 -*-
"""Setup file used to compile EponaSolver.pyx
Created on Mon Jun 04 11:16:29 2012

@author: Christopher Strickland
"""
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np

ext_modules = [Extension("EponaGSolverC",["EponaGSolverC.pyx"],\
include_dirs=[np.get_include()])]

setup(
    name = 'Epona Solver module',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)

#to build, run "python setup.py build_ext --inplace" in a terminal
#or get C code + type annotation by running cython -a EponaSolverC.pyx
#
#With SDK tools, you will first need to set the environment. Run:
#"C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\bin\vcvars64.bat"
#set MSSdk=1
#set DISTUTILS_USE_SDK=1
#
#Note: replace with vcvars32.bat for 32-bit version
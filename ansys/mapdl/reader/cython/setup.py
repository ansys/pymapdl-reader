"""Installation file for ansys-mapdl-reader-support"""
from io import open as io_open
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext
import os
import platform
import re
import struct
import subprocess
import sys

import numpy as np

# Facilities to install properly on Mac using clang
def is_clang(bin):
    proc = subprocess.Popen([bin, '-v'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    output = str(b'\n'.join([stdout, stderr]).decode('ascii', 'ignore'))
    return not re.search(r'clang', output) is None


class build_ext(_build_ext):
    """ build class that includes numpy directory """

    def build_extensions(self):
        if os.name != 'nt':
            binary = self.compiler.compiler[0]
            if is_clang(binary):
                for e in self.extensions:
                    e.extra_compile_args.append('-stdlib=libc++')

                    if platform.system() == 'Darwin':
                        # get the minor version
                        mac_version, _, _ = platform.mac_ver()
                        minor = [int(n) for n in mac_version.split('.')][1]

                        # libstdc++ is deprecated in recent versions of XCode
                        if minor >= 9:
                            e.extra_compile_args.append('-mmacosx-version-min=10.9')
                            e.extra_compile_args.append('-stdlib=libc++')
                            e.extra_link_args.append('-mmacosx-version-min=10.9')
                            e.extra_link_args.append('-stdlib=libc++')
                        else:
                            e.extra_compile_args.append('-mmacosx-version-min=10.7')
                            e.extra_link_args.append('-mmacosx-version-min=10.7')

        _build_ext.build_extensions(self)


def compiler_name():
    """ Check compiler and assign compile arguments accordingly """
    import re
    import distutils.ccompiler
    comp = distutils.ccompiler.get_default_compiler()
    getnext = False

    for a in sys.argv[2:]:
        if getnext:
            comp = a
            getnext = False
            continue
        # separated by space
        if a == '--compiler' or re.search('^-[a-z]*c$', a):
            getnext = True
            continue
        # without space
        m = re.search('^--compiler=(.+)', a)
        if m is None:
            m = re.search('^-[a-z]*c(.+)', a)
            if m:
                comp = m.group(1)

    return comp


# Assign arguments based on compiler
compiler = compiler_name()
if compiler == 'unix':
    cmp_arg = ['-O3', '-w']
else:
    cmp_arg = ['/Ox', '-w']


# Get version from version info
__version__ = "0.0.dev0"
install_requires = [
    'numpy>=1.16.0'
]

# perform python version checking
# this is necessary to avoid the new pip package checking as vtk does
# not support Python 32-bit as of 17 June 2021.
is64 = struct.calcsize("P") * 8 == 64
if not is64:
    raise RuntimeError('\n\n``ansys-mapdl-reader`` requires 64-bit Python.\n'
                       'Please check the version of Python installed at\n'
                       '%s' % sys.executable)

setup(
    name='ansys-mapdl-reader-support',
    version=__version__,
    description='Support library for ansys.mapdl.reader',
    long_description=open('README.rst').read(),
    long_description_content_type='text/x-rst',
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'License :: OSI Approved :: MIT License',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Operating System :: MacOS',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    url='https://github.com/pyansys/pymapdl-reader/ansys/mapdl/reader/cython',

    # Build cython modules
    cmdclass={'build_ext': build_ext},
    include_dirs=[np.get_include()],
    ext_modules=[
                 Extension('ansys.mapdl.reader._archive',
                           ['_archive.pyx',
                            'archive.c'],
                           extra_compile_args=cmp_arg,
                           language='c',),

                 Extension('ansys.mapdl.reader._reader',
                           ['_reader.pyx',
                            'reader.c',
                            'vtk_support.c'],
                           extra_compile_args=cmp_arg,
                           language='c',),

                 Extension("ansys.mapdl.reader._relaxmidside",
                           ["_relaxmidside.pyx"],
                           extra_compile_args=cmp_arg,
                           language='c'),

                 Extension("ansys.mapdl.reader._cellqual",
                           ["_cellqual.pyx"],
                           extra_compile_args=cmp_arg,
                           language='c'),

                 Extension("ansys.mapdl.reader._binary_reader",
                           ["_binary_reader.pyx",
                            "binary_reader.cpp"],
                           extra_compile_args=cmp_arg,
                           language='c++'),
                 ],

    python_requires='>=3.7.*',
    keywords='vtk MAPDL ANSYS cdb full rst',
    install_requires=install_requires,
)

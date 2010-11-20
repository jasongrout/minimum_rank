from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

# IMPORTANT: this package *requires* Sage, so should be installed using sage -python
# See the Makefile for an example command
import os
SAGE_ROOT=os.environ["SAGE_ROOT"]

include_dirs=[
    SAGE_ROOT+'/local/include/csage/',
    SAGE_ROOT+'/local/include/',
    SAGE_ROOT+'/devel/sage/',
    ]

setup(
    # I'm not sure exactly what the following should be
    # I should probably look this up in some distutils guide or something
    #name = "minrank",
    #description
    #version = "1.0",
    #author = "Jason Grout and others",
    #author_email = "jason.grout@drake.edu",
    #url = "github...";,
    #description = "Minimum Rank",
    #packages=["minimum_rank"],
    #package_dir={"minimum_rank":"src"},

    # classifiers=[
    #   'Development Status :: 4 - Beta',
    #   'Environment :: X11 Applications :: GTK',
    #   'Intended Audience :: End Users/Desktop',
    #   'Intended Audience :: Developers',
    #   'License :: OSI Approved :: GNU General Public License (GPL)',
    #   'Operating System :: POSIX :: Linux',
    #   'Programming Language :: Python',
    #   'Topic :: Desktop Environment',
    #   'Topic :: Text Processing :: Fonts'
    #   ]

    ext_modules = [
        Extension(
            "zero_forcing_wavefront", # name of extension
            ["zero_forcing_wavefront.pyx"], # filename of Cython source
            include_dirs=include_dirs,
            ),

        Extension(
            "zero_forcing_64", # name of extension
            ["zero_forcing_64.pyx"], # filename of Cython source
            include_dirs=include_dirs,
            ),

        Extension(
            "Zq_c", # name of extension
            ["Zq_c.pyx"], # filename of Cython source
            include_dirs=include_dirs,
            ),


        # Extra options that could be specified in an extension tuple
        #language="c++",              # this causes Cython to create C++ source
        #libraries=["stdc++", ...],             # ditto
        #extra_link_args=[...],       # if needed
        
        ],
    
    # Standard stuff
    cmdclass = {'build_ext': build_ext},
)

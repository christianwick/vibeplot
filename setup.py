"""
Usage:
    python setup.py sdist
    python setup.py py2exe
"""

import os
import glob
from distutils.core import setup
import platform
import matplotlib

if platform.system() == "Windows":
    import py2exe

OPTIONS = dict(
    includes="matplotlib numpy sip openbabel".split(),
    excludes="""numpy.random scipy wx
                _imagingtk PIL._imageingtk
                ImageTk PIL.ImageTk FixTk
                Tkconstants Tkinker tcl""".split(),
    dll_excludes="MSVCP80.dll MSVCR80.dll MSVCP90.dll".split(),
    skip_archive=False,
    #compressed=2,
    bundle_files=3,)

if platform.system() == "Windows":
    data_files = matplotlib.get_py2exe_datafiles()
    data_files.append(
        (".",
         glob.glob(os.path.join(os.environ["BABEL_DATADIR"], "../*.obf"))))
        # may also need {libstdinchi,libxml2,iconv}.dll
else:
    data_files = None

setup(
    name="QVibePlot",
    version="0.14.0",
    author="Mathias Laurin",
    author_email="Mathias.Laurin+vibeplot@gmail.com",
    url="http://vibeplot.sf.net",
    license="Python Software Foundation License",
    packages="vibeplot vibeplot.utils".split(),
    #package_dir={'vibeplot': 'src'},
    #package_data=package_data,
    data_files=data_files,
    requires="maplotlib (>=1.1);numpy;sip;PyQt4;openbabel".split(";"),
    scripts=["QVibePlot.py"],
    windows=[dict(script="QVibePlot.py")],
    zipfile=None,
    options=dict(py2exe=OPTIONS),
    classifiers=os.linesep.join(
        s for s in """
        Development Status :: 4 - Beta
        License :: OSI Approved :: Python Software Foundation License
        Intended Audience :: Education
        Intended Audience :: Science/Research
        Environment :: X11 Applications :: Qt
        Operating System :: POSIX :: Linux
        Operating System :: Microsoft :: Windows :: Windows NT/2000
        Operating System :: POSIX :: BSD
        Programming Language :: Python :: 2
        Topic :: Education
        Topic :: Scientific/Engineering :: Chemistry
        Topic :: Scientific/Engineering :: Physics
        Topic :: Scientific/Engineering :: Visualization
        """.splitlines() if s)
)

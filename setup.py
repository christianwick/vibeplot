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
    datadir = os.environ["BABEL_DATADIR"]
    data_files.append(  # may be empty if openbabel was compiled from source
        (".",
         glob.glob(os.sep.join((datadir, "..", "*.obf")))))
    data_files.append(
        (".",  # supporting dlls for openbabel
         [os.sep.join((datadir, "..", dll_name + ".dll"))
          for dll_name in "iconv libinchi libxml2 xdr-0".split()]))
else:
    data_files = []

data_files.append(
    (os.sep.join("data orbimol".split()),
     glob.glob(os.sep.join("data orbimol *.freq".split()))))
data_files.append(
    ("licenses", glob.glob(os.sep.join("licenses *.txt".split()))))

setup(
    name="QVibePlot",
    version="0.14.2",
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

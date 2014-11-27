"""
Usage:
    python setup.py sdist
    python setup.py py2exe
"""

from __future__ import print_function
import os
import sys
import glob
from distutils.core import setup
import platform
import matplotlib
from vibeplot import __version__

if platform.system() == "Windows":
    import py2exe

OPTIONS = dict(
    includes="matplotlib numpy sip openbabel".split(),
    excludes="""scipy wx
                _imagingtk PIL._imageingtk
                ImageTk PIL.ImageTk FixTk
                Tkconstants Tkinker tcl""".split(),
    dll_excludes="MSVCP80.dll MSVCR80.dll MSVCP90.dll".split(),
    skip_archive=False,
    optimize=2,
    compressed=2,
    bundle_files=2,)


def isfile(path_):
    _isfile = os.path.isfile(path_)
    print(path_, _isfile, file=sys.stderr, sep="\t")
    return _isfile


if platform.system() == "Windows":
    # matplotlib
    data_files = matplotlib.get_py2exe_datafiles()

    # Open Babel
    ob_libdir = os.environ.get("BABEL_LIBDIR", "../ob-build/bin")  # plugins
    plugins = [
        os.path.join(ob_libdir, plugin + ".obf")
        for plugin in  # selected plugins for openbabel
        """
        APIInterface
        acesformat crystal09format gamessformat gamessukformat
        gaussformat moldenformat molproformat mopacformat
        nwchemformat orcaformat qchemformat vaspformat cmlformat
        plugin_ops
        """.split()]
    dlls = [
        os.path.join(ob_libdir, dll + ".dll")
        for dll in  # supporting dlls for openbabel
        "iconv libinchi libxml2 xdr-o".split()]
    data_files.append((".", [plugin for plugin in plugins if isfile(plugin)]))
    data_files.append((".", [dll for dll in dlls if isfile(dll)]))

    # PyQt5
    import PyQt5
    data_files.append((".",
                       [os.path.join(next(PyQt5.__path__.__iter__()),
                                     "libEGL.dll")]))
    import PyQt5.plugins.platforms as platforms
    data_files.append(("platforms",
                       [os.path.join(next(platforms.__path__.__iter__()),
                                     "qwindows.dll")]))
else:
    data_files = []

data_files.append((os.path.join("data", "orbimol"),
                   glob.glob(os.path.join("data", "orbimol", "*.freq"))))
data_files.append(
    ("licenses", glob.glob(os.sep.join("licenses *.txt".split()))))

setup(
    name="QVibePlot",
    version=__version__,
    author="Mathias Laurin",
    author_email="Mathias.Laurin+vibeplot@gmail.com",
    url="http://vibeplot.sf.net",
    license="Python Software Foundation License",
    packages=". vibeplot vibeplot.utils".split(),
    #package_dir={'vibeplot': 'src'},
    #package_data=package_data,
    data_files=data_files,
    requires="maplotlib (>=1.4);numpy;sip;PyQt4;openbabel".split(";"),
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

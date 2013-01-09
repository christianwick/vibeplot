"""
Usage:
    python setup.py sdist
    python setup.py py2exe
"""

from distutils.core import setup
import platform
import matplotlib

if platform.system() == "Windows":
    import py2exe

OPTIONS = dict(includes=['matplotlib', 'numpy', 'sip'],
               excludes=['scipy',
                         'pybel', 'openbabel',
                         'Tkconstants', 'Tkinter', 'tcl',
                         'numpy.random' ],
               dll_excludes=['MSVCP90.dll'],
               skip_archive=False,
               compressed=2,
               bundle_files=1,)

data_files = (matplotlib.get_py2exe_datafiles() 
              if platform.system() == "Windows" else None)
#package_data = {'': ['vibeplot/UI/*']}

setup(
    name="QVibePlot",
    version="0.13.1",
    author="Mathias Laurin",
    author_email="Mathias.Laurin+vibeplot@gmail.com",
    url="http://vibeplot.sf.net",
    license="Python Software Foundation License",
    packages=['vibeplot', 'vibeplot.utils', 'vibeplot.sdg',
              'vibeplot.parser', 'oasa', 'oasa.graph', ],

    #package_dir={'vibeplot': 'src'},
    #package_data=package_data,
    data_files=data_files,
    requires=["matplotlib (>=1.1)", "numpy", "sip", "PyQt4"],
    provides=["oasa (0.13.1)"],
    scripts=["QVibePlot.py"],
    windows=[dict(script="QVibePlot.py")],
    zipfile = None,
    options=dict(py2exe=OPTIONS),
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: Python Software Foundation License',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Environment :: X11 Applications :: Qt',
        'Operating System :: POSIX :: Linux',
        'Operating System :: Microsoft :: Windows :: Windows NT/2000',
        'Operating System :: POSIX :: BSD',
        'Programming Language :: Python :: 2',
        'Topic :: Education',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
)

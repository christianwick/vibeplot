vibeplot application programming interface (API)
================================================

Interactive session with vibeplot
---------------------------------

The vibeplot library may be used directly to generate the plots interactively
with Matplotlib or in batch.  The example provided here is for the 0.14 series.

>>> import pybel
>>> import openbabel as ob
>>> import matplotlib.pyplot as plt
>>> from vibeplot.plotter import VibrationPlotter
>>> import vibeplot.sdg as sdg

We first instanciate a :class:`VibrationPlotter` object and provide is with a
figure and axes.

>>> vb = VibrationPlotter()
>>> vb.fig = plt.figure()
>>> vb.axes = vb.fig.add_subplot(111)

We now load the data from a file named ``benzene.molf`` in the Molden format
[#molden.fmt]_.

>>> molecule = pybel.readfile("data/benzene.input")
>>> vb.molecule = molecule.OBMol
>>> vibrations = ob.toVibrationData(molecule.OBMol.GetData(ob.VibrationData))
>>> vb.normal_coordinates = vibrations.GetLx()

The 'VibrationPlotter' object created above contains a series of
'VibrationPlotter.get_*' methods that are called in order to draw a molecule
and add the markers.  Extra keyword arguments passed to them is forwarded to
Matplotlib (an example could be to set the 'linewidth' or 'zorder').  Moreover,
two convenience methods are provided, 'plot_molecule' and 'plot_vibration',
that call the relevant 'get_*' with good default arguments.

>>> vb.plot_molecule()

We now have the skeleton of a benzene molecule. We pop and display a random
vibration.

>>> vb.plot_vibration(5)
>>> plt.show()

And let us save the results using another convenience method 'save_molecule'.

>>> vb.save_molecule("%.0f.png" % frequency)

This simple demonstration shows most of what is possible with vibeplot.

.. image:: 3100.png
   :scale: 30 %

.. [#molden.fmt] Description of the `Molden format`_
.. [#mpl.col] API of `Matplotlib.collections`_

.. _Molden format: http://www.cmbi.ru.nl/molden/molden_format.html
.. _Matplotlib.collections: http://matplotlib.sourceforge.net/api/collections_api.html



.. automodule:: vibeplot.plotter

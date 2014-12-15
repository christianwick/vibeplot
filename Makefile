doc:
	sphinx-build -b html sphinx html

doc-all:
	sphinx-build -aE -b html sphinx html

upload-html:
	rsync -avzP -e ssh html/ mathias_laurin@web.sourceforge.net:/home/project-web/vibeplot/htdocs/ 

profile:
	python -m cProfile -s cumulative QVibePlot.py

plop:
	python -m plop.collector QVibePlot.py

lint:
	pylint --include-ids=y --disable=C0103 QVibePlot.py

rcc:
	pyuic5 -o qvibeplot_ui.py qvibeplot.ui
	pyrcc4 -o rcc4.py rcc.qrc
	pyrcc5 -o rcc5.py rcc.qrc

distribution:
	python setup.py sdist

exe:
	python setup.py py2exe

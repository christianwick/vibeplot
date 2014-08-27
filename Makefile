doc:
	sphinx-build -b html sphinx html

doc-all:
	sphinx-build -aE -b html sphinx html

upload-html:
	rsync -avzP -e ssh html/ mathias_laurin@web.sourceforge.net:/home/project-web/vibeplot/htdocs/ 

profile:
	python -m cProfile -s cumulative QVibePlot.py

lint:
	pylint --include-ids=y --disable=C0103 QVibePlot.py

distribution:
	python setup.py sdist

exe:
	python setup.py py2exe

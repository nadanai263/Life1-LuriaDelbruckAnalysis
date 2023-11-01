# Life1-LuriaDelbruckAnalysis
Scripts for analysis of yeast fluctuation experiments. 

<img src="/plots/class_histLDWhites.png" width="480" alt="LD distribution" align="center" hspace="40" vspace="15">

Part of the first year Life 1 course at the School of Biological Sciences, University of Edinburgh. Feedback and corrections to Nadanai Laohakunakorn ([nadanai.laohakunakorn@ed.ac.uk](mailto:nadanai.laohakunakorn@ed.ac.uk)).

### Instructions
There are many ways to run this. The easiest is to copy the Jupyter notebook (`LDAnalysis.ipynb`) to [Google Colab](https://colab.research.google.com/), along with the `data` folder. Make sure you have the file `sampledata.xlsx` present inside the `data` folder. 

If you would like to run locally on your computer, clone the repository to your hard drive, and make sure you have a working Python installation. Install the packages from the pipfile (ideally use an environment manager such as [conda](https://anaconda.org/) or [pipenv](https://pipenv.pypa.io/en/latest/)). Then start Jupyter notebook by using the command line to navigate to the repository, and running 

	jupyter notebook

You can also run the analysis directly on the command line,  by running

	python LDAnalysis.py

The Jupyter noteook is preferred as it is more interactive.

Finally, you can also run through a [Docker](https://www.docker.com/) container. With Docker installed on your system, navigate to your repository and start the Jupyter notebook by typing

	docker run -p 8888:8888 --rm -it -v "$PWD":/home/jovyan jupyter/datascience-notebook

The first time you do this it will take some time as Docker installs the required container.

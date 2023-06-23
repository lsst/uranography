Installation
============


Quick Installation
------------------

Appearing shortly: ``pip install uranography``

coming soon: ``conda install -c conda-forge uranography``


For Developer Use
-----------------

First, clone the `uranography` repository:

::

 git clone git@github.com:lsst/uranography.git
 cd uranography


Create a conda environment for it:

::

 conda create --channel conda-forge --name uranography --file requirements.txt python=3.11


If you want to run tests, install the test requirements as well:

::

 conda activate uranography
 conda install -c conda-forge --file=test-requirements.txt


Install the `uranography` project into this environment (from the uranography directory):

::

 pip install -e .


In order to make uranography available in jupyter,
install the kernel from the new jupyter environment:

::

 python -m ipykernel install --user --name=uranography

If you need to remove the environment and/or kernel:

::

  conda env remove -n uranography
  jupyter kernelspec uninstall uranography


Building Documentation
----------------------

An online copy of the documentation is available at https://uranography.lsst.io,
however building a local copy can be done as follows:

::

 pip install "documenteer[guide]"
 cd docs
 package-docs build


The root of the local documentation will then be `docs/_build/html/index.html`.

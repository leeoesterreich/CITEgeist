Installation
============

Prerequisites
-------------

CITEgeist requires **Python 3.10 or higher** and the following dependencies:

* `numpy <https://numpy.org/>`_
* `pandas <https://pandas.pydata.org/>`_
* `scanpy <https://scanpy.readthedocs.io/>`_
* `squidpy <https://squidpy.readthedocs.io/>`_
* `scipy <https://scipy.org/>`_
* `anndata <https://anndata.readthedocs.io/>`_
* `scikit-learn <https://scikit-learn.org/>`_
* `matplotlib <https://matplotlib.org/>`_
* `networkx <https://networkx.org/>`_
* `tqdm <https://tqdm.github.io/>`_
* `joblib <https://joblib.readthedocs.io/>`_
* `esda <https://pysal.org/esda/>`_
* `libpysal <https://pysal.org/libpysal/>`_

GPU Requirement (cuOPT)
-----------------------

CITEgeist's cell-type proportion solver (Module 3) uses **quadratic programming (QP)** via the
`cuOPT <https://developer.nvidia.com/cuopt-logistics-optimization>`_ library, which requires:

* An **NVIDIA GPU** with CUDA support (8 GB+ VRAM recommended)
* The `cuOPT` Python package — **not available on PyPI**

Install cuOPT via the NVIDIA package index or a pre-built RAPIDS/NGC container::

   pip install cuopt-cu12 --extra-index-url https://pypi.nvidia.com

There is **no CPU fallback** for the proportion solver. CPU-only environments cannot run
Module 3.  All other modules (GEX deconvolution, program discovery, integration) work
without a GPU.

Installation
------------

CITEgeist is published on PyPI::

   pip install citegeist

To install from source instead::

   git clone https://github.com/leeoesterreich/CITEgeist.git
   cd CITEgeist
   pip install -e .

.. note::

   The package name is ``citegeist`` (lowercase), but the importable package
   name is ``CITEgeist`` (capital C). After installation::

      import CITEgeist
      from CITEgeist import CitegeistModel

For a development installation with test and lint dependencies:

.. code-block:: bash

   pip install -e .[dev]

Verification
------------

To verify your installation:

.. code-block:: python

   import CITEgeist
   print(f"CITEgeist {CITEgeist.__version__} installed successfully!")

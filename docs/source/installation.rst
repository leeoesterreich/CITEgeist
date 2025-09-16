Installation
============

Prerequisites
-------------

CITEgeist requires Python 3.7 or higher and the following dependencies:

* `numpy <https://numpy.org/>`_
* `pandas <https://pandas.pydata.org/>`_
* `scanpy <https://scanpy.readthedocs.io/>`_
* `scipy <https://scipy.org/>`_
* `pyarrow <https://arrow.apache.org/docs/python/>`_
* `gurobipy <https://www.gurobi.com/>`_ (requires Gurobi license)

Installation
------------

Install CITEgeist from source:

.. code-block:: bash

   git clone https://github.com/leeoesterreich/CITEgeist.git
   cd CITEgeist
   pip install -e .

Gurobi Setup
------------

CITEgeist uses Gurobi for optimization. You'll need:

1. A valid Gurobi license
2. Gurobi installed on your system
3. The `gurobipy` Python package

For detailed Gurobi installation instructions, see the `Gurobi documentation <https://www.gurobi.com/documentation/>`_.

Verification
------------

To verify your installation:

.. code-block:: python

   import citegeist_model
   print("CITEgeist installed successfully!")

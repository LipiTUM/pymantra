.. MANTRA's documentation master file, created by
   sphinx-quickstart on Thu Apr  7 16:23:34 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. |PyPI| |License| |Documentation| |Codecov| |Downloads|

pymantra - Reaction-centered Metabolic Network Analysis
=======================================================================

.. image:: _static/Figure1.svg
   :width: 80%
   :align: center

**mantra** is a conceptual approach to compute estimates for the change of
metabolic reaction activity between two groups of samples. It relies on linear
relationships between substrate and product metabolites of a reaction and how
the coefficients of these relationships change between conditions. In addition
to analyzing metabolomics data, mantra also provides a correlation-based
approach for *multi-omics integration*.

As an approach to provide smaller, mechanistically interpretable results based
on both the reaction estimates (and multi-omics associations) and the metabolic
network structure, network enrichment on the basis of a simulated-annealing
assisted local search is used.

The ``pymantra`` package provides all functionalities for computing changes in
reaction activity, multi-omics associations and performing the network
enrichment as well as reporting and plotting their results. Additionally, it
contains utilities to perform metabolite ID mapping.

The general workflow of the package is summarized below.

.. mermaid:: /workflow.mmd


Manuscript
----------
If you would like to learn more about the details of the methodology and see
some real-world results please check out our manuscript :cite:`koehler2023`

Getting Started
===============

Installation
------------

To get instructions on installation, follow our
:doc:`installation guide <installation>`.

Examples
--------

For an introduction to mantra's functionality and how to use it, check out our
:doc:`examples <examples/index>`.

API Documentation
-----------------

The full documentation is found in the :doc:`API section <api/index>`

.. toctree::
   :maxdepth: 2
   :caption: Contents

   installation
   examples/index
   api/index

.. toctree::
   :maxdepth: 2
   :caption: About

   contact
   citation
   license
   GitHub <https://github.com/lipitum/pymantra>


References
==========
.. bibliography::
    :cited:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. |PyPI| image:: https://img.shields.io/pypi/v/pymantra.svg
    :target: https://pypi.org/project/pymantra
    :alt: PyPI

.. |License| image:: https://img.shields.io/pypi/l/pymantra.svg
    :target: https://pypi.org/project/pymantra
    :alt: License

.. |Downloads| image:: https://pepy.tech/badge/pymantra
    :target: https://pepy.tech/project/pymantra
    :alt: Downloads

.. |Docs|  image:: https://img.shields.io/readthedocs/pymantra
    :target: https://pymantra.readthedocs.io/en/stable
    :alt: Documentation

.. |Codecov| image:: https://codecov.io/gh/pymantra/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/pymantra
    :alt: Coverage

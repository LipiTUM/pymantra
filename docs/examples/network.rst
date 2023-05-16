.. _network_generation:

Network Generation
==================

..
    TODO: include remote API call

We start by defining a set of metabolites and microbial species, for which we
want to extract the metabolic network.

Please note that all names need to be given as internal names/IDs. For more
information on how to map them see :ref:`name_mapping`.

.. literalinclude::
    ../../examples/network.py
    :language: python
    :lines: 1-13
    :lineno-start: 1

To generate a valid network from a local neo4j database, we use the
``NetworkGenerator`` class. This requires neo4j to be installed (we recommend
using docker for this) and the mantra.dump file to be moved in the correct
place (see the installation manual for more details). To initialize the
``NetworkGenerator`` object you only need to pass the URI and (optionally)
a user and password for authentication.

.. literalinclude::
    ../../examples/network.py
    :language: python
    :lines: 15
    :lineno-start: 15

Subsequently, we extract all edges for the given metabolites and microbial
species from the database. The edges are returned by edge type. All reactions
extracted with the given `reaction_organism` are either human-catalysed
reactions or catalysed by one of the microbial species passed.

.. literalinclude::
    ../../examples/network.py
    :language: python
    :lines: 17-19
    :lineno-start: 17

Finally, we can use the extracted edges to generate a ``networkx.DiGraph``
object, which automatically contains the correct node- and edge-type annotation
required for downstream analyses.

Additionaly, we call ``reduce_reaction_nodes`` on the resulting graph to avoid
having multiple reaction nodes with the same substrates and products, as there
is no way to distinguish them in the downstream computations.

.. literalinclude::
    ../../examples/network.py
    :language: python
    :lines: 21-23
    :lineno-start: 21

Full Example Code
~~~~~~~~~~~~~~~~~
.. literalinclude:: ../../examples/network.py
    :language: python
    :lineno-start: 1

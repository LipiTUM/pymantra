Database Module
===============

The database module handles the generation of metabolic networks
from a given list of metabolites and optionally genes and/or microbial
species.

Depending on whether the database is accessed through the provided REST
framework or directly through neo4j (regardless of whether you are using neo4j
in a docker container or outside) different query classes should be used.

For the REST API, please use ``APINetworkGenerator``, for neo4j queries
``NetworkGenerator``.

Online/API Database Class
~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: pymantra.database.APINetworkGenerator
    :members:

Local Database Class
~~~~~~~~~~~~~~~~~~~~
.. autoclass:: pymantra.database.NetworkGenerator
    :members:


Functions
~~~~~~~~~
.. autofunction:: pymantra.database.reduce_reaction_nodes

Exceptions
~~~~~~~~~~

.. autoclass:: pymantra.database.IncorrectEdgeType
.. autoclass:: pymantra.database.IncorrectNodeType

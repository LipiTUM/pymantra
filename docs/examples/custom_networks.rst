.. _custom_network:

Using custom networks
=====================

In addition to using the mantra database as a basis for enrichment analysis,
it is also possible to use self-generated networks. The only restriction for
these networks is that they need to have the same node- and edge-types
annotated with attributes named 'node_type' and 'edge_type'.

We recommend to also set the additional attributes for compatibility with the
provided plotting functions:

* reaction: Formula, Description
* metabolite: Name
* organism: nodeLabel
* gene: nodeLabel


The following names are currently used:

* nodes
    * reaction
    * metabolite
    * organism
    * gene
* edges
    * SUBSTRATE (metabolite - reaction)
    * PRODUCT (reaction - metabolite)
    * REACTION_ORGANISM (reaction - organism)
    * REACTION_GENE (reaction - gene)

We start with a set of metabolite and reaction nodes (and optionally organism
and gene nodes) and edges.

.. literalinclude::
    ../../examples/custom_network.py
    :language: python
    :lines: 1-21
    :lineno-start: 1

The first step is to add the nodes with the required 'node_type' attribute.

.. literalinclude::
    ../../examples/custom_network.py
    :language: python
    :lines: 23-30
    :lineno-start: 23

Subsequently, the edges can be added and the edge type can directly be inferred
using the static variable `EDGE_BY_NODE_TYPE`, which is a dictionary, where the
keys are 2-tuples of source and target node types.

To avoid having multiple reactions with the same substrates and products (these
will always get the same reaction activity change), we use the
``reduce_reaction_nodes`` function.

.. literalinclude::
    ../../examples/custom_network.py
    :language: python
    :lines: 32-38
    :lineno-start: 32

Finally, we can plot the generated graph. This is also a partial sanity check,
if node colours and shapes need to match the match the ones in the legend.

.. literalinclude::
    ../../examples/custom_network.py
    :language: python
    :lines: 31-32
    :lineno-start: 31

The graph can then be used to compute reaction models and multi-omics
associations and to run a local search as in :ref:`reaction_enrichment`
and :ref:`metabolome_microbiome_enrichment`. The only requirement is of
course that the metabolite, organsim and gene names match between graph and
data.

Full Example Code
~~~~~~~~~~~~~~~~~

.. literalinclude::
    ../../examples/custom_network.py
    :language: python
    :lineno-start: 1
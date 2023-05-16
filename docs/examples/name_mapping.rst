.. _name_mapping:

Name Mapping
============

When using mantra's neo4j database to generate a metabolic network matching
experimental data, you will need to convert feature names or their database
IDs to internal IDs.

Microbial Organisms
~~~~~~~~~~~~~~~~~~~

Microbes are simply accessed by their species name, e.g.
"Bacteroides uniformis". Importantly, make sure names start with a capital
letter and that separation is happening through whitespaces and not underscores
or similar.

Metabolites
~~~~~~~~~~~

Recommended database IDs are KEGG and HMDB, followed by Reactome and Virtual
Metabolic Human. Some other database like ChEBI or NCBI are supported, but are
likely to give less matches.

For mapping metabolite names and database IDs, we provide some functions.
Their usage is shown in the following code snippets.

We start by loading the required packages/functions and load a pre-defined set
of metabolites, specified by their "common" name.

.. literalinclude::
    ../../examples/name_mapping.py
    :language: python
    :lines: 1-11
    :lineno-start: 1

Next, we query the database IDs for these common names. For this, mantra uses
the Metaboanalyst API. In case you already have database IDs for your
metabolites, you can skip this step.

.. literalinclude::
    ../../examples/name_mapping.py
    :language: python
    :lines: 13-16
    :lineno-start: 13

Lastly, we use an internal database to convert from (in this case) HMDB IDs to
mantra IDs. For all available database sources options please see the
documentation of `NameMapper`.

.. literalinclude::
    ../../examples/name_mapping.py
    :language: python
    :lines: 18-24
    :lineno-start: 18


Full Example Code
~~~~~~~~~~~~~~~~~

.. literalinclude::
    ../../examples/name_mapping.py
    :language: python
    :lineno-start: 1
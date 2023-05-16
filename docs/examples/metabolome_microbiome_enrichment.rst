.. _metabolome_microbiome_enrichment:

Reaction-based Metabolome-Microbiome Integration
================================================

We first load some example data from github. This is data used also used in
the manuscript. We simply load a .csv file with metabolite measurements and
a group annotation file, together with a pre-computed metabolite-reaction
graph. For more details on how to generate such a graph and what its structure,
looks like see :ref:`network_generation` and :ref:`custom_network`. The
metabolite data should contain the samples in rows.

.. literalinclude:: ../../examples/metabolite_microbiome_enrichment.py
    :language: python
    :lines: 1-14
    :lineno-start: 1

Once all data is loaded, we first compute the linear models and the difference
in residuals between samples groups by calling ``compute_reaction_estimates``
and add them to the graph in the required way using ``add_reaction_estimates``.
For how to use custom reaction models, see :ref:`custom_reaction_model`.

.. literalinclude:: ../../examples/metabolite_microbiome_enrichment.py
    :language: python
    :lines: 17-19
    :lineno-start: 17

Next, the associations between the reaction values and the microbiome
measurements are computed by ``compute_microbiome_associations`` and add them
to the graph by ``add_microbiome_associations``. The correlations and p-values
can also be analyzed seperately.

.. literalinclude:: ../../examples/metabolite_microbiome_enrichment.py
    :language: python
    :lines: 21-24
    :lineno-start: 21

Finally, we can initialize a ``MultiOmicsLocalSearch`` object and run the
local search. Afterwards, the local search object contains a number of
functions to report and plot the results.

.. literalinclude:: ../../examples/metabolite_microbiome_enrichment.py
    :language: python
    :lines: 26-44
    :lineno-start: 26


Entire Example Code
~~~~~~~~~~~~~~~~~~~
.. literalinclude:: ../../examples/metabolite_microbiome_enrichment.py
    :language: python
    :lineno-start: 1

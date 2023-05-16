.. _reaction_enrichment:

Finding dysregulated metabolic reactions
========================================

We first load some example data from the packages.  We simply load a .csv file
with metabolite measurements and a group annotation file, together with a
pre-computed metabolite-reaction graph. For more details on how to generate
such a graph and what its structure, looks like see :ref:`network_generation`
and :ref:`custom_network`. The metabolite data should contain the samples in
rows.

.. literalinclude:: ../../examples/metabolite_enrichment.py
    :language: python
    :lines: 1-12
    :lineno-start: 1

Once all data is loaded, we compute the linear models and the difference in
residuals between samples groups by calling ``compute_reaction_estimates`` and
add them to the graph in the required way using ``add_reaction_estimates``.
For how to use custom reaction models, see :ref:`custom_reaction_model`.
For including confounding variables seek :ref:`confounder_correction`.

.. literalinclude:: ../../examples/metabolite_enrichment.py
    :language: python
    :lines: 14-16
    :lineno-start: 14

Finally, we can initialize a ``MetaboliteLocalSearch`` object and run the
local search. Afterwards, the local search object contains a number of
functions to report and plot the results.

.. literalinclude:: ../../examples/metabolite_enrichment.py
    :language: python
    :lines: 18-26
    :lineno-start: 18


Full Example Code
~~~~~~~~~~~~~~~~~
.. literalinclude:: ../../examples/metabolite_enrichment.py
    :language: python
    :lineno-start: 1

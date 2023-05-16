.. _confounder_correction:

Including confounding variables in dysregulation analysis
=========================================================

We first load some example data from the packages.  We simply load a .csv file
with metabolite measurements and a group annotation file, together with a
pre-computed metabolite-reaction graph. For more details on how to generate
such a graph and what its structure, looks like see :ref:`network_generation`
and :ref:`custom_network`. The metabolite data should contain the samples in
rows. In addition to the example data, we also generate some confounding data.

.. literalinclude:: ../../examples/estimation_with_correction.py
    :language: python
    :lines: 1-18
    :lineno-start: 1

Once all data is loaded, we compute the linear models and the difference in
residuals between samples groups by calling ``compute_reaction_estimates`` and
add them to the graph in the required way using ``add_reaction_estimates``.
The confounder correction is happening inside ``compute_reaction_estimates``.
In this case we are using the following formula:
``products ~ substrates + site + gender + (1|age)``

.. literalinclude:: ../../examples/estimation_with_correction.py
    :language: python
    :lines: 20-25
    :lineno-start: 20

After this step the enrichment can be performed the same way as demonstrated in
:ref:`reaction_enrichment`. Furthermore, the multi-omics association function
``compute_multiomics_associations`` (see :ref:`metabolome_microbiome_enrichment`
for reference) also allows for correction in the same way as
``compute_reaction_estimates``.

Full Example Code
~~~~~~~~~~~~~~~~~
.. literalinclude:: ../../examples/estimation_with_correction.py
    :language: python
    :lineno-start: 1

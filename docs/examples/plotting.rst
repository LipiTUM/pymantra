.. _plotting:

Plotting mantra results
=======================

The plotting example uses mantra's standard example dummy dataset, which is
loaded first. Additionally, a custom color palette is defined to showcase how
to set custom colors.

.. literalinclude:: ../../examples/plotting.py
    :language: python
    :lines: 1-18
    :lineno-start: 1

First the full example graph is plotted using ``plot_directed_graph``. A
function for undirected plotting with the same interface is also available.

.. literalinclude:: ../../examples/plotting.py
    :language: python
    :lines: 21-25
    :lineno-start: 21

Next the plotting function for estimated reaction values is conducted. The
created plot is a violin plot indicating the residual distributions of the
compared sample groups

.. literalinclude:: ../../examples/plotting.py
    :language: python
    :lines: 27-32
    :lineno-start: 27

Last, the plotting of multi-omics associations is done. After computing the
associations ``plot_correlation_differences`` is used to plot a heatmap showing
the differences in associations for each reaction/multi-omic feature pair.
In this example ``set_zero`` is set to ``False``, due to the example data not
having many significant associations, but in practice we recommend using the
default ``True``.

Additionally, the correlations between reaction/multi-omic feature pairs with
the highest absolute difference in correlations between groups can be plotted
as individual scatter plots to get a more detail view on their changes.

.. literalinclude:: ../../examples/plotting.py
    :language: python
    :lines: 34-51
    :lineno-start: 47


Full Example Code
~~~~~~~~~~~~~~~~~
.. literalinclude:: ../../examples/plotting.py
    :language: python
    :lineno-start: 1

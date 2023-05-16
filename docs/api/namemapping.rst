Name Mapping
============

The name mapping module implements all functionalities required to map external
IDs like HMDB, ChEBI, KEGG, InCHI keys etc. to the internal IDs used in the
mantra database.

While it is possible to use the base classes for mapping between certain
databases, we recommend using the wrapper class (`NameMapper`).

Wrapper
~~~~~~~
.. autoclass:: pymantra.namemapping.NameMapper
    :members:


Base Classes
~~~~~~~~~~~~
.. autoclass:: pymantra.namemapping.ChEBIQuery
    :members:
.. autoclass:: pymantra.namemapping.HMDBQuery
    :members:
.. autoclass:: pymantra.namemapping.ReactomeQuery
    :members:
.. autoclass:: pymantra.namemapping.MantraDBQuery
    :members:


Exceptions
~~~~~~~~~~
.. autoclass:: pymantra.namemapping.UnknownMappingError

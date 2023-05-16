Installation
============

Since installing ``pymantra`` requires compiling C++ extensions, we recommend
installing with pip, which uses binary (pre-compiled) wheels. These are
available for most Linux, MacOS and Windows architectures.

Installation using pip
~~~~~~~~~~~~~~~~~~~~~~

To install from ``pip`` run one of the following options ::

    # install with package dependencies
    pip install pymantra
    # for including the dependencies to run the paper experiments
    pip install pymantra[experiments]
    # for including the development dependencies
    pip install pymantra[dev]
    # for including the documentation dependencies
    pip install pymantra[docs]

Installation from source
~~~~~~~~~~~~~~~~~~~~~~~~

To install from source first clone the github repository including submodules ::

    git clone https://github.com/lipitum/pymantra.git --recursive


Make sure you have a C++ compiler installed (recommended are gcc for
Linux and OS X and Visual C++ for Windows). Please make sure the
compiler matches the one, with which your python distribution was installed.

The only c++ dependency for compilation is
the `boost library <https://robots.uc3m.es/installation-guides/install-boost.html>`_.
In case it is already installed on your system (with version >= 1.77) you can
also drop the ``--recursive`` flag in the ``git clone`` call. In addition to
installing boost from source you can also use ``conda``.

If you want to use mantra's multiprocessing option, you need to have
OpenMP version > 2 available on your system. On macOS this might require to
install a different compiler or using homebrew to install the library
(``brew install libomp``). On Windows the only option to use OpenMP (at the
time of writing) is using the Windows Subsystem Linux (WSL).

Windows user need to have the Visual Studio Code C++ development kit installed
see e.g. `here <https://stackoverflow.com/questions/60322655/how-to-use-cython-on-windows-10-with-python-3-8>`_
for help.

To compile and install navigate to the mantra folder and call ::

    pip install .

Similar to the regular ``pip`` installation you can also install the optional
dependencies.

To verify your installation run all unit tests ::

    python -m pytest

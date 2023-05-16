# Contribution Guidelines

## Future of the Project

There are a few things that are planned to be implemented in the future

* [ ] multivariate non-parametric test to compute significances for
	  multi-variate residual comparison
* [ ] correcting reaction estimates by substrate concentrations of
	  'incoming' reactions
* [ ] Linear-mixed model based multi-omics associations
* [ ] gene ID mapping

Of course any additional features are wellcome! If you have any ideas or
requests, please either open up an issue or contact us via e-mail first.
The same of course applies to question about the current TODOs.
We're happy discussing anything you are interested in!

## Minimal Contribution Recommendations

All newly implemented functions should try to stay as close as possible to
the existing API. Modification of existing functions should aim at being
backwards compatible.

Submission of contributions should be in the form of pull requests,
pass code convention checks, have documentation and have unit-tests
(see below), if they implemented new features or functions.

Before submitting, please make sure both html and pdf generation of sphinx
documentation are working.

## Conventions

### Documentation

Most importantly, document your code well. To keep documentation in a
consistent manner and to be able to easily generate documentation homepages,
please stick to the [NumPy docstring format](
https://numpydoc.readthedocs.io/en/latest/format.html).

In addition to documentation, which should help those using the package,
please also comment your code, unless it's self-explanatory (that of course
requires well-named variables etc.).


### Error Messages

Please be mindful that users of this package might not be very experienced
with python. Therefore, in addition to simple API design, try to write clear
and meaningful error messages, that do not necessarily require users to fully
understand your code.

Also, not all no-brainers are no-brainers. Just keep [this meme](
https://twitter.com/code_memez/status/1239789862091292680) in mind.


### Code Style
In order to keep code style in a somewhat consistent manner, we stick to the
[PEP 8 conventions](https://www.python.org/dev/peps/pep-0008/). To ensure
these, please check your code (e.g. before commiting) using
[`flake8`](https://flake8.pycqa.org/).

```bash
flake8 test.py
```

In many IDEs, you can also set these conventions automatically.


### Unit testing

Please make sure all your main functions have unit tests.
Ideally, these should contain manually generated 'gold-standards' (when
possible) to ensure correct functionality.
If this is not possible, they should at least guarantee that all main functions
run through without errors and that edge cases are handled properly.

For automated testing, we use [pytest](https://docs.pytest.org/en/7.2.x/).
Simply running

```bash
python -m pytest
```

in the base folder of this repository will automatically find run all tests
in the 'test' folder, as long as their file names start or end with 'test'.
After installing the package the command can be shortened to `pytest`.
It is also possible to run individual files. For most IDEs there also exists
a pytest integration.


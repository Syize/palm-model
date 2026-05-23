# Developer guide

Python development environment, formatting, code structure and editor configuration

---

This document provides a brief overview of the code, its style and the automated testing of `palm_csd`.

First, it is recommended to create and activate a virtual environment in the root folder of `palm_csd` repository, either manually with

```bash
python -m venv .venv
source .venv/bin/activate
```

or using editor-specific tools. In order to install all the necessary dependencies for development, run the following command in the root directory of the repository:

```bash
pip install -r requirements-dev.txt
```

The set-up of these tools is done in the `pyproject.toml` file. The checks and tests detailed below run for each Merge Request on PALM's Gitlab system as set-up in the `.gitlab-ci.yml` file. All checks and tests have to succeed before new code can be merged.

## Formatting and checks

### Python code formatting, linting and type checks

The code is formatted and linted using [`Ruff`](https://docs.astral.sh/ruff/) to ensure a consistent code style. Furthermore, the code should include [type hints](https://docs.python.org/3/library/typing.html), which are checked using [`mypy`](https://mypy-lang.org/). Currently, we do not impose any further style guidelines. In order to format the code, run the following command in the root directory of the repository:

```bash
ruff format
```

or use

```bash
ruff format --check
```

to check if the code is formatted correctly. Note that the sorting of the imports is handled by the linter functionality described next.

We also use the code linter functionality of `ruff` to check for potential issues in the code. To run the linter, use the following command:

```bash
ruff check
```

In order to sort the imports automatically, run

```bash
ruff check --fix --select "I"
```

In order to fix all fixable issues, run

```bash
ruff check --fix
```

Make sure to check the result of this.

We use the Google style ([here](https://google.github.io/styleguide/pyguide.html#s3.8-comments-and-docstrings) and [here](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html)) for docstrings. They are also checked with the above command.

The static type checker `mypy` can be run with the following command:

```bash
mypy .
```

### Functionality tests

The functionality of `palm_csd` is tested using the `pytest` framework. The tests are located in the `tests` directory. To run the tests, execute the following command in the root directory of the repository:

```bash
pytest
```

The test coverage is measured using the `pytest-cov` extension and printed after running the tests.

### Markdown linting

The `README.md` and documentation is linted using [pymarkdown](https://github.com/jackdewinter/pymarkdown). To lint the markdown files, run the following command in the root directory of the repository:

```bash
pymarkdown scan -r README.md docs/
```

## Code structure and functionality

The main executable is `palm_csd.py` in the main folder of the repository. It calls `create_driver` in `palm_csd/create_driver.py`. This routine initializes various objects, reads in the user configuration, and handles the debug option. It further creates for each domain a `CSDDomain` object. The main calculations happen in a loop over all domains. Here, several other functions are called for each domain.

Variables and functions that start with an underscore `_` are considered private and are intended to be only used within their module. Variables that are written in all capital letters are considered constants.

Members of `Enum` classes are deliberately written in lowercase. The names are used for user input and output in the netCDF files. Using uppercase would result in a more complicated code that has to use `lower()` and `upper()` several times.

### Multi-dimensional arrays

The main data structure used for multi-dimensional fields are [Numpy](https://numpy.org/)'s [`MaskedArray`](https://numpy.org/doc/stable/reference/maskedarray.html)s. A `MaskedArray` stores both, the values in the [`data`](https://numpy.org/doc/stable/reference/maskedarray.baseclass.html#numpy.ma.MaskedArray.data) and a mask in the [`mask`](https://numpy.org/doc/stable/reference/maskedarray.baseclass.html#numpy.ma.MaskedArray.mask) attribute. The mask is a boolean array usually with the same shape as the data array. If an element in the mask is `True`, the corresponding value is considered invalid. The advantage of this approach is that invalid values can be masked out for any [data type (`dtype`)](https://numpy.org/doc/stable/reference/arrays.dtypes.html), while the [`numpy.nan`](https://numpy.org/doc/stable/reference/constants.html#numpy.nan) value can only be used with float data. Note that if all elements of the mask are `False`, the mask might be represented by [`nomask`](https://numpy.org/doc/stable/reference/maskedarray.baseclass.html#numpy.ma.nomask), which corresponds to a single `np.False_` value. `np.False_` works similarly as a Python `False` value. If a mask with the same shape as the data array is required, use [`getmaskarray`](https://numpy.org/doc/stable/reference/generated/numpy.ma.getmaskarray.html) instead of directly accessing the `mask` attribute. This is _not_ required as arguments for [`ma.where`](https://numpy.org/doc/stable/reference/generated/numpy.ma.where.html) because its arguments are broadcasted as needed. Also `.all()` and `.any()` methods as well as `~`, `&` and `|` operators work on single `np.False_` values.

### Configuration and validation

The yaml user configuration is read in by the `__init__` of the `CSDConfig` class in `palm_csd/csd_config.py`. The different sections of the configuration are dealt with by the different `CSDConfig*` classes inheriting from `CSDConfigElement`. This class, in turn, inherits from [pydantic](https://docs.pydantic.dev/)'s [`BaseModel`](https://docs.pydantic.dev/latest/concepts/models/).

pydantic offers extensive validation of the input data in terms of data type and data ranges. For example, instance variables declared as `float` will be `float`, even if an `int` value is supplied; a variable declared as a `Path` will convert a `str` to a `Path` object. With declaring a variable as [`Field(default=..., ge=..., le=...)`](https://docs.pydantic.dev/latest/concepts/fields/), a default, a minimum and a maximum value can be set. Furthermore, custom validators can be defined. A [`@field_validator` function](https://docs.pydantic.dev/latest/concepts/validators/#field-validators) can be used to check the validity of one variable, while a [`@model_validator` function](https://docs.pydantic.dev/latest/concepts/validators/#model-validators) can be used to check the validity of the whole model. These custom validators can be applied with `mode="before"` to check the raw input data before pydantic's standard validation. `mode="after"` validators are applied after pydantic's standard validation. The validators can be applied with the decorators given above. Alternatively, custom types can be declared with [`Annotated`](https://docs.pydantic.dev/latest/concepts/validators/#annotated-validators) together with `BeforeValidator` and `AfterValidator`. The latter approach has the advantage that the custom type can be applied to different variables.

The default, minimum and maximum values are defined in `palm_csd/data/value_defaults.csv`. This file is read by `_populate_defaults`.

In `CSDConfigElement`, a counter how often an instance is created is stored. Depending on `_unique`, more than one instance allowed, which corresponds to several configuration sections of the same type. For example, while only one `setting` section is allowed in the yaml configuration, multiple `domain` sections are allowed. The `model_validator` `_validate_unique_counter` checks that.

### Logging

`palm_csd` uses the [Python `logging` module](https://docs.python.org/3/library/logging.html) for message output, which is mainly set-up in `palm_csd/__init__.py`. While we use simple terminal output, the code uses extensively the different available log levels. These levels allow selective output according to the selected output level. In addition to the default levels, we introduce the level `STATUS`. With this, the levels, in that order, `DEBUG`, `STATUS`,  `INFO`, `WARNING`, `ERROR` and `CRITICAL` are available. `STATUS` is used for messages that indicate what part of the processing is currently done. `INFO` gives more detailed information to the processing step. `WARNING` is used for messages that indicate potential problems or automatic data modifications. `ERROR` is currently used and `CRITICAL` is meant for errors that stop the processing. `DEBUG` messages give additional messages that are not normally not needed. In order to support the `STATUS` level, the `Logger` class is extended with the `status` method in `StatusLogger`.

In order to emphasize the different log levels, different colors and indentions are used. The colors and indention levels are applied in the `ColorFormatter`. By default, `STATUS` messages are not indented and the other messages are indented by one level. If a custom indention level is required, the `StatusLogger` offers the `*_indent()` methods. This is mainly used to collect messages together: If a new message gives additional information to the message before, the new message is indented by one level relative to the former message.

The root logger is initialized in `palm_csd/__init__.py` and the default log level is set to `STATUS` with `logger_root.setLevel(STATUS)`. This means that only messages with the level `STATUS` and higher are displayed. In each module, a module logger is created with `logger = logging.getLogger(__name__)`. The log level can be adjusted for each module individually. In particular, depending on the user setting, the respective log level is set to `DEBUG` with `logger.setLevel(logging.DEBUG)` in `create_driver` in `palm_csd/create_driver.py`.

All logging messages end with a `.` if they include at least a verb and a noun. This applies to virtually all messages.

### Object-oriented programming

The code is written in an object-oriented manner where deemed useful. In particular, this approach is used to collect variables and functions that belong together in a class. For example, for each domain, a `CSDDomain` object is created. This object contains all required dimensions (of type `NCDFDimension`) and variables (of type `NCDFVariable`) of this domain. The different configuration sections are represented by different classes as described [above](#configuration-and-validation). Each single tree is represented by a `DomainTree` object. There are also objects that focus on computational routines and only store their required parameters. Examples are `GeoConverter` and `CanopyGenerator` objects.

### Memory management

The main program code in `create_driver` is divided into several functions. This is done to enhance the readability but also to reduce the memory footprint of the program. Python frees memory of objects that are not accessible anymore using its garbage collector. Every object that is available only in a function is freed after the function is left. Thus, in most routines, the result fields are not directly stored in the `CSDDomain` domain object and its `NCDFVariable`'s `values` attribute but in local variables. The content of these variables is stored in the output NetCDF before the function is left.

## Editor configuration

### Visual Studio Code

[Visual Studio Code](https://code.visualstudio.com/) supports most of the tools used in this project. With the [Python extension](https://marketplace.visualstudio.com/items?itemName=ms-python.python) and the [Ruff extension](https://marketplace.visualstudio.com/items?itemName=charliermarsh.ruff), you can format the code and run the linter directly in the editor. In the `settings.json`, set

```json
"[python]": {
    "editor.defaultFormatter": "charliermarsh.ruff",
    "editor.formatOnSave": true,
}
```

to format the code with Ruff on save. The `settings.json` file can be accessed by pressing `Ctrl + Shift + P` and selecting the `Preferences: Open User Settings (JSON)`.

The imports can be sorted automatically by the Ruff extension by pressing `Ctrl + Shift + P` and `Organize Imports`. This can also be done on save by adding

```json
    "editor.codeActionsOnSave": {
        "source.organizeImports": "explicit"
    }
```

to the `[python]` settings above.

A virtual environment can be created by pressing `Ctrl + Shift + P` and `Python: Create Environment...`. Select the appropriate environment type (if in doupt `Venv`), Python interpreter and select `requirements-dev.txt` as the requirements file. This should create a virtual environment, install all required dependencies and activate the environment.

`palm_csd`'s automatic tests can be run from the Test Explorer. Just select the pytest framework here. In order to debug the Python code, install the [Python Debugger extension](https://marketplace.visualstudio.com/items?itemName=ms-python.debugpy).

Currently, there is no extension available for `pymarkdown`. There is, however, [an extension for `markdownlint`](https://marketplace.visualstudio.com/items?itemName=DavidAnson.vscode-markdownlint) that features similar rules as `pymarkdown`.

With the [Gitlab extension](https://marketplace.visualstudio.com/items?itemName=GitLab.gitlab-workflow), you can create and manage Merge Requests directly from Visual Studio Code. The extension also supports Gitlab's CI/CD pipeline.

# ProPhyle – Developers

## Contributing Guidelines

### ProPhyle subcommands

* All steps of downstream analysis that are expected to be done by the end users should be implemented through subcommands (e.g., `prophyle index`). Please, keep the main script, `prophyle`, as simple as possible.
* All other (auxiliary) procedures should be available only through prophyle\_.py scripts.


### CLI

* Only short command-line parameters should be used (e.g., `-l`, not `--log-file`).
* CLI parameters with a capital letter are switchers without argument (e.g., `-R`).
* When possible, all required program arguments should be passed through positional arguments.


### Source codes

* Each program or script should contain the name of the author(s) and the license (ideally MIT).
* Python scripts should follow the [Google Python Style Guide](https://google.github.io/styleguide/pyguide.html).


### Testing

* A-tests (small) are run on Travis. B-tests (big) are run only locally (they are too long to be tested after each commit).
* Every program should (ideally) have a unit test.


### Packaged Python scripts

* Scripts used by prophyle should be named prophyle\_\*
* Auxiliary testing scripts should not be included in the ProPhyle package.

## Releasing

* Every release info should be structured as "New - Improvements - Fixes"
* Releases should be available from PyPI, BioConda and Github

## FAQs

> Travis tests don't pass due to missing packages even though everything seems to be fine.

Try to remove Travis caches (see the button _More options_).


# ProPhyle – Developers

## Contributing Guidelines

### ProPhyle subcommands

* The steps which are supposed to be done by the end user should be implemented through subcommands.
* All other procedures should be implemented through prophyle\_.py scripts.


### CLI

* Only short command-line parameters should be used (e.g., `-L`).
* CLI parameters with a capital letter are switchers without arguments (e.g., `-R`).
* When possible, required program arguments should be passed through positional arguments


### Source codes

* Each program or script should contain a name of the author and the license (ideally MIT).
* Python scripts should follow the [Google Python Style Guide](https://google.github.io/styleguide/pyguide.html).


### Testing

* A-tests (small) are run on Travis. B-tests (big) are run only locally (they are too long to be tested after each commit).
* Every program should (ideally) have a unit test.


### Packaged Python scripts

* Scripts used by prophyle should be named prophyle\_\*
* Auxiliary scripts (e.g., those used only in the tests) should not be included in the ProPhyle package


## FAQs

> Travis tests don't pass due to missing packages even though everything seems to be fine.

Try to remove Travis caches (see the button _More options_).

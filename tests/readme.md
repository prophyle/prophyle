# Tests

* A - unit tests
* B - small scale integration tests
* C - large scale integration tests

## General patterns

* temporary files should start by _, repo files shouldn't
* diff files should start with __
* diff should be called with the -c param
* sam outputs should pass through samtools view -h, which checks basic syntax
* make clean should clean all non-repo files

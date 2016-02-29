
# CONTRIBUTING #

We encourage users, researchers and programmers to contribute to this
free and open source project. This can be achieved by reporting bugs and
commiting code to the github repository
https://github.com/yhoogstrate/segmentation-fold. To streamline and
archive communication in an univocal way, we encourge contributors to
only use this channel, Github, to contribute to segmentation-fold.
To contribute to segmentation-fold, please make use the following
documentation standard:
[http://www.stack.nl/~dimitri/doxygen/index.html](http://www.stack.nl/~dimitri/doxygen/index.html).

## Tidy ##

All source-code is formatted using Asteric Style: [http://astyle.sourceforge.net/](http://astyle.sourceforge.net/).
The corresponding configuration file is available at '[share/.astylerc](https://raw.githubusercontent.com/yhoogstrate/segmentation-fold/master/share/.astylerc)'.
If you want to contribute to the code you should have it installed to 'beautify' the code, by simply running:

	make tidy

Cmake will make use of the correct syntax file and beautify all C++ and hpp files in
[src/](https://github.com/yhoogstrate/segmentation-fold/tree/master/src),
[include/](https://github.com/yhoogstrate/segmentation-fold/tree/master/include),
and
[test/](https://github.com/yhoogstrate/segmentation-fold/tree/master/test).

## Tests ##

To simplify the reviewing process of submitted code the project contains
unit and functional tests. These tests have to be passed in order to get
a positive review. These tests also inspect memory leaks using valgrind.
To run the test on (your copy of the) code before doing a pull request, run:

	cmake .
	make clean
	make readme
	make
	make test
	ctest -V -T memcheck

This will re-build the readme, re-compile the code, and does testing with and
without memory leak checking. If you can't get it working but you still believe
your change is worth submitting, don't worry. Whenever you do a pull request
TravisCI will automatically run the tests for you.

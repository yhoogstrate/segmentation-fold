
# INSTALLATION #

## Prerequisites ##

Segmentation-fold does **not** install out of the box because it depends
on two additional libraries and a tool required for building:

	cmake
	boost library (-system)
	boost library (-test)
	boost library (-xml)

In Debian and Ubuntu you can install these dependencies with the following command:

	sudo apt-get install cmake libboost-system-dev libboost-test-dev libboost-filesystem-dev

In Arch linux you can install these packages with the folllowing command:

	sudo pacman -S git cmake boost boost-libs

## Recommended packages ##

To create the corresponding documentation, you should have installed the
following package:

	doxygen (>= 1.8.3)

The doxygen package is version specific because Markdown support is
implemented in 1.8.3 and above.

To make the style of the syntax of (novel) code automatically identical
to existing code, you should have installed the package:

	astyle

You can install the recommended packages in Ubuntu or Debian with:

	sudo apt-get install doxygen astyle

You can install the recommanded packages in Arch with:

	sudo pacman -S doxygen astyle

Details on running doxygen and astyle are given in the sections
[Get the documentation](https://github.com/yhoogstrate/segmentation-fold#get-the-documentation)
and
[CONTRIBUTING](https://github.com/yhoogstrate/segmentation-fold#contributing).

## Build and install ##

After having the prerequisites installed, you can downadload, compile
and install segmentation-fold with following commands:

	 $ git clone https://github.com/yhoogstrate/segmentation-fold
	 $ cd segmentation-fold
	 $ cmake -DCMAKE_BUILD_TYPE=release .
	 $ make
	 $ make test
	 $ sudo make install

In case you do not have administrator rights on a particular machine,
you can install as follows:

	 $ git clone https://github.com/yhoogstrate/segmentation-fold
	 $ cd segmentation-fold
	 $ cmake -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=~/.local .
	 $ make
	 $ make test
	 $ make install

You will then find the binary in the directory:
	<home directory>/.local/bin/

If you want to run segmentation-fold for anything else than development,
you should compile with the -DCMAKE_BUILD_TYPE=release flag because by
default compiling is done with the 'debug' flag. The two differences
between debug and release are that (1) debug contains certain
(unneccesairy if all code behaves as expected) checks on variables and
(2) has no optimized, but faster, compilation. Therefore, debug mode
will produce a binary that has a slower performance, but a higher level
of safety.

If you are using a compiled version of segmentation-fold and you don't
known the original compilation settings, you should run:

	segmentation-fold --version

This will report the version number, the git commit-SHA1 and the build
type.

## Get the documentation ##

Run the following commands to get the latest version of the
documentation:

	cmake .
	make readme
	make doc

The first command, "cmake ." creates the doxygen description file from
its template and sets e.g. the package version variable. Then "make
readme" creates the latests version of README.md based on the Markdown
files present in this project. Then "make doc" creates the documentation
from all the doxygen comments within the code, including the latest
README.md as main page.

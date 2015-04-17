
# INSTALLATION #

## Prerequisites ##

Segmentation-fold does **not** install on most systems of itself because
it depends on two additional libraries and an installation library.

	cmake
	boost library (-test)
	boost library (-xml)

In ubuntu you can install them very easily with the following command:

	sudo apt-get install libboost-test-dev libboost-filesystem-dev

## Recommended packages ##

To create the corresponding documentation, you should have the following
package installed:

	doxygen (>= 1.8.3)

The doxygen package is version specific because of the Markdown
integration. If you want to change code you should have the following
tool installed to 'beautify' the code:

	astyle

## Build and install ##

Please follow the instructions given in the '[INSTALL](https://github.com/yhoogstrate/segmentation-fold/blob/master/INSTALL)' file, e.g. by
using one of the following commands:

	 $ git clone https://github.com/yhoogstrate/segmentation-fold
	 $ cd segmentation-fold
	 $ cmake -DCMAKE_BUILD_TYPE=release .
	 $ make
	 $ make test
	 $ sudo make install

In case you do not have administrator rights on a particular machine, you can
install as follows:

	 $ git clone https://github.com/yhoogstrate/segmentation-fold
	 $ cd segmentation-fold
	 $ cmake -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=~/.local .
	 $ make
	 $ make test
	 $ make install

You will then find the binary in the directory:
	<home directory>/.local/bin/

## Get the documentation ##

Run the following commands to get the latest version of the
documentation:

	cmake .
	make readme
	make doc

The first command, "cmake ." creates the doxygen description file from
its template and sets e.g. the package version variable. Then
"make readme" creates the latests version of README.md based on the Markdown
files present in this project. Then "make doc" creates the documentation
from all the doxygen comments within the code, including the latest
README.md as main page.

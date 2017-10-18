#!/bin/bash

# This file is still a work in progress, and is very UConn specific at the moment!

module load gcc/4.8.5
module load boost/1.61.0
module load eigen/3.3.2

# remove any existing distr directory and create new empty one
rm -rf distr
mkdir distr

# create variable to store top-level directory
cd distr
DISTRDIR="$( pwd )"

# create populate distr with source files
cd $DISTRDIR
cp ../arbor/arbor/*.hpp .
cp ../arbor/arbor/*.cpp .

# create Makefile
echo "arbor : main.o arbor.o" > Makefile
echo "	g++ -L$BOOST_LIBRARYDIR -o arbor main.o arbor.o -lm -lboost_program_options -lboost_regex" >> Makefile
echo "arbor.o : arbor.cpp arbor.hpp grapedb.hpp grape.hpp lot.hpp utilfunc.hpp xarbor.hpp" >> Makefile
echo "	g++ -std=c++11 -O2 -I$BOOST_INCLUDEDIR -I$EIGEN_INCLUDEDIR -c arbor.cpp" >> Makefile
echo "main.o : main.cpp arbor.hpp grapedb.hpp grape.hpp lot.hpp utilfunc.hpp xarbor.hpp" >> Makefile
echo "	g++ -std=c++11 -O2 -I$BOOST_INCLUDEDIR -I$EIGEN_INCLUDEDIR -c main.cpp" >> Makefile
echo "clean :" >> Makefile
echo "	rm arbor main.o arbor.o" >> Makefile

# run make
make

# create bash script to run arbor
echo "#!/bin/bash" > grape-arbor.sh
echo "module load boost/1.61.0" >> grape-arbor.sh
echo "module load gcc/4.8.5" >> grape-arbor.sh
echo "arbor $@" >> grape-arbor.sh
chmod +x grape-arbor.sh

echo "done!"


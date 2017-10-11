#!/bin/bash

# This file is still a work in progress, and is very UConn specific at the moment!

# The Boost C++ Libraries were successfully built!
#
# The following directory should be added to compiler include paths:
#
#     /home/plewis/software/boost_1_65_1
#
# The following directory should be added to linker library paths:
#
#     /home/plewis/software/boost_1_65_1/stage/lib

module load gcc/4.8.5
module load boost/1.61.0
module load eigen/3.3.2
module list
g++ --version

#export BOOST_INCLUDEDIR="/home/plewis/software/boost_1_65_1"
#export BOOST_LIBRARYDIR="/home/plewis/software/boost_1_65_1/stage/lib"
#export EIGEN_INCLUDEDIR="/commmon/opt/eigen/3.3.2"
export EIGEN_INCLUDEDIR="/home/plewis/software"

# remove any existing distr directory and create new empty one
rm -rf distr
mkdir distr

# create variable to store top-level directory
cd distr
DISTRDIR="$( pwd )"

# create Makefile.am inside distr
cd $DISTRDIR
echo "AUTOMAKE_OPTIONS = foreign" > Makefile.am
echo "SUBDIRS = src" >> Makefile.am

# create distr/src directory and populate with source files
cd $DISTRDIR
mkdir src
cd src
cp ../../arbor/arbor/*.hpp .
cp ../../arbor/arbor/*.cpp .
echo "CXX = g++" > Makefile.am
echo "AM_CPPFLAGS = -Wall -O2 -std=c++11 -I$EIGEN_INCLUDEDIR -I$BOOST_INCLUDEDIR -L$BOOST_LIBRARYDIR" >> Makefile.am
echo "LIBS = -lboost_program_options -lboost_regex" >> Makefile.am
echo "bin_PROGRAMS = arbor" >> Makefile.am
echo "arbor_SOURCES = arbor.cpp main.cpp" >> Makefile.am

# run autoscan to create default configure.ac
cd $DISTRDIR
autoscan
mv configure.scan configure.ac

# use sed to specify program name, version, and contact email
cd $DISTRDIR
sed -i 's/\[FULL-PACKAGE-NAME\]/arbor/' configure.ac
sed -i 's/\[VERSION\]/1.0/' configure.ac
sed -i 's/\[BUG-REPORT-ADDRESS\]/paul.lewis@uconn.edu/' configure.ac
#sed -i 's|AC_OUTPUT|AC_OUTPUT(Makefile src/Makefile)|' configure.ac
sed -i '/AC_INIT/ a\
AM_INIT_AUTOMAKE(arbor, 1.0)' configure.ac

echo "running aclocal..."

# run aclocal to generate aclocal.m4 containing macros for automake (e.g. AM_INIT_AUTOMAKE)
aclocal

echo "running autoheader..."

# run autoheader to generate config.h.in needed by automake
autoheader

echo "running automake..."

# run automake to produce Makefile.in for every Makefile.ac
# The argument --add-missing provides default scripts for reporting errors, installing, etc.
automake --add-missing

echo "running autoconf..."

# run autoconf to create the configure script
autoconf

# run configure
cd $DISTRDIR
./configure
make

echo "done!"


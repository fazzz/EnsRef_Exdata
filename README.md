# EnsRef_Exdata

## Requirement

libLBFGS (http://www.chokkan.org/software/liblbfgs/)

## Installation

``
./configure CPPFLAGS=-ILBFGS-INSTALL-DIR/include LDFLAGS=-LLBFGS-INSTALL-DIR/lib
make
make install
``

## How to use

``
ensResEx experimentaldata simulationdata weigh optimalobservalvalues
``

## Options

[--Del] `error width`

[--Expm] `parameter file name`

[--h] `show the 


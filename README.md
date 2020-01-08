# EnsRef_Exdata

## Requirement

libLBFGS (http://www.chokkan.org/software/liblbfgs/)

## Installation

``
./configure CPPFLAGS=-I*LBFGS-INSTALL-DIR*/include LDFLAGS=-L*LBFGS-INSTALL-DIR*/lib
make
make install
``

## How to use

``
ensResEx *experimental-data-name* *simulation-data-name* *weights-output-name* *optimal-observal-values-output-name*
``

## Options

[--Del] `error width`

[--Expm] `parameter file name`

[--h] `show the help messeage`




# EnsRef_Exdata

Tool to re-weight the snapshots of molecular simulation to fit to the experimental data

## Requirement

libLBFGS (http://www.chokkan.org/software/liblbfgs/)

## Installation

``
./configure CPPFLAGS=-ILBFGS-INSTALL-DIR/include LDFLAGS=-LLBFGS-INSTALL-DIR/lib
make && make install
``

## How to use

``
ensResEx *experimental-data-name* *simulation-data-name* *weights-output-name* *optimal-observal-values-output-name*
``

## Options

[--Del] `error width`

[--Expm] `parameter file name`

[--h] `show the help messeage`




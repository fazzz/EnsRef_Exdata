# EnsRef_Exdata

This is a tool to re-weight the snapshots of molecular simulation to fit to the experimental data.

## The algorithms

In this method, the weighing factors of each snapshot of MD simulation are reweighted
for the improvement of the agreement between experimental measurements and
the corresponding calculated values based on simulation data without undesired without
undesired changes of the distribution generated by the simulation. For the aim,
the optimal weighting factors are derived
by the maximizing the following objective function L,

<img src="https://latex.codecogs.com/gif.latex?L_{n}(\{&space;\omega&space;\},\rho_{n})&space;=D&space;(\rho&space;\mid&space;\rho^{\mathrm&space;sim.})&space;&plus;&space;\lambda&space;(\sum_{i=1}^{N}&space;\omega_{i}&space;-&space;1)&space;&plus;&space;\rho_{n}&space;([&space;\max&space;(0,\chi^{2}(\omega)&space;-&space;\delta^{2})&space;])^{2}" />

where <img src="https://latex.codecogs.com/gif.latex?\inline&space;\rho" />
and <img src="https://latex.codecogs.com/gif.latex?\inline&space;\rho^{sim.}" /> 
are the optimized and original distribution,
<img src="https://latex.codecogs.com/gif.latex?\inline&space;D(\rho&space;\mid&space;\rho^{sim.})" />
is the Kullback-Leibler divergence
between <img src="https://latex.codecogs.com/gif.latex?\inline&space;\rho" />
and <img src="https://latex.codecogs.com/gif.latex?\inline&space;\rho^{sim.}" />, 
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\lambda" />
is Lagrange multiplier,
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\chi^2\left(\omega\right)" />
is the square error of the given experimental measurement and calculated correspondence from
simulation trajectory whose snapshots are 
weighted by factors 
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\left\{\omega_i\right\}" />,
N is the number of the snapshots, <img src="https://latex.codecogs.com/gif.latex?\inline&space;\rho_{n}" /> is the penalty parameter for enforcing <img src="https://latex.codecogs.com/gif.latex?\inline&space;\chi^2\left(\omega\right)" /> is
smaller than some given tolerance parameter <img src="https://latex.codecogs.com/gif.latex?\inline&space;\delta^2" />.

The desired optimal weighting factors <img src="https://latex.codecogs.com/gif.latex?\inline&space;\left\{\omega_i\right\}" />
that satisfy the condition <img src="https://latex.codecogs.com/gif.latex?\inline&space;\chi^2\left(\omega\right)\le\delta^2" />
is obtained after sequential solving of
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\left\{\mathcal{L}_1(\rho_1),\cdots,\mathcal{L}_M(\rho_M)\right\}" />, 
where <img src="https://latex.codecogs.com/gif.latex?\inline&space;\rho_1<\cdots<\rho_M" />.

## Requirement

libLBFGS (http://www.chokkan.org/software/liblbfgs/) is required for the optimization.

## Installation

For installation of `ensRefEx` it is required to install `libLBFGS` correctly.
We assume that you have installed `libLBFGS` in the directory named $LBFG*.

Firstly, in the directry where `Makefile.in` exsist, type as following

``
./configure CPPFLAGS=-I$LBFG/include LDFLAGS=-$LBFG/lib
make && make install
``

The generated program name will be `ensRefEx`.

## How to use

If you have experimental data and simulation data,
you can run `ensRefEx` as follows:

``
ensResEx *experimental-data* *simulation-data* *weights-output* *optimal-observal-values-output*
``

where `experimental-data` is the file of experimental data, and `simulation-data` is the file of simulation data,
and `weights-output` and `optimal-observal-values-output` are the files of output data.

The examples of each file can be seen in `examples`.

## Options

The followings are the options for `ensRefEx`.

[--Del] `error width` (default)

[--Expm] `parameter file name` (default)

[--h] `show the help messeage`


## Parameter files

The parameter file describes parameter r for exterior-point method.
The example is as follows:

``
1.0e-2 1.0e-1 1.0e0 1.0e1 1.0e2
``

## Examples

This repository contains an example of how `ensRefEx` works.
The first example of the method is a lysozyme simulation
(PDB ID: 6LYT ) in explicit water molecules.
In the case , the MD trajectory was reweighted according to
the experimental values of J3-coupling.
J3-coupling, the spin-spin interaction of 1H and another 1H via H-C-H bond, is one of the
most important NMR parameters for the determination of protein structure
in solution. The dependence of J3-coupling constants on protein dihedralangles are
described by Karplus equation.

To run the examples completley, you must install `python3` and `mdtraj`.

To generate the experimental values vs. compputed values W/WO reweighing:

``
cd examples;
./Fig2.sh
``

To generate the relationship between <img src="https://latex.codecogs.com/gif.latex?\inline&space;D(\rho&space;\mid&space;\rho^{sim.})" /> 
and <img src="https://latex.codecogs.com/gif.latex?\inline&space;\chi^{2}" />
with multiple <img src="https://latex.codecogs.com/gif.latex?\inline&space;\delta^{2}" />.

``
./Fig3.sh
``

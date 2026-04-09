# ReaclibFit

ROOT-based tool for fitting thermonuclear reaction rates in the REACLIB 7-parameter form.

## Build

c++ reaclib_fit.cpp -o reaclib_fit $(root-config --cflags --libs)

## Run

./reaclib_fit -c exp.conf

./reaclib_fit -c exp.conf --input ExpRate_central.dat --output exp_central.fit --no-gui

## Input Format

The fitter reads the first two columns only:

T9   reaction_rate

Extra columns after the second one are ignored.

## Directories

conf/: config files and term definitions
input/: input data and reference fit files
output/: fit output files
macros/: ROOT macros

## Resonance Convention

For resonant terms, the code uses:

Er = Ex - Threshold

## Macros

Generate central/lower/upper tables:
root 'macros/make_rate_variants.C("input/ExpRate.dat")'

Plot multiple fit results:
root 'macros/plot_fit_collection.C("input/ExpRate_central.dat,input/ExpRate_lower.dat,input/ExpRate_upper.dat","output/exp_central.fit,output/exp_lower.fit,output/exp_upper.fit")'

## Reference

Cyburt et al., ApJS 189 (2010) 240


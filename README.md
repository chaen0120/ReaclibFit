# ReaclibFit

ROOT-based tool for fitting thermonuclear reaction rates in the REACLIB 7-parameter form.

## Build

- Compile with:
  ```bash
  c++ reaclib_fit.cpp -o reaclib_fit $(root-config --cflags --libs)
  ```

## Run

- Run with GUI:
  ```bash
  ./reaclib_fit -c exp.conf
  ```

- Run without GUI:
  ```bash
  ./reaclib_fit -c exp.conf --input ExpRate_central.dat --output exp_central.fit --no-gui
  ```

## Input Format

- The fitter reads the first two columns only:
  ```text
  T9   reaction_rate
  ```

- Extra columns after the second one are ignored.

## Directories

- `conf/`: config files and term definitions
- `input/`: input data and reference fit files
- `output/`: fit output files
- `macros/`: ROOT macros

## Resonant Terms

- For resonant terms, the code uses:
  ```text
  Er = Ex - Threshold
  ```

## Macros

- Generate central/lower/upper tables:
  ```bash
  root 'macros/make_rate_variants.C("input/ExpRate.dat")'
  ```

- Plot multiple fit results together:
  ```bash
  root 'macros/plot_fit_collection.C("input/ExpRate_central.dat,input/ExpRate_lower.dat,input/ExpRate_upper.dat","output/exp_central.fit,output/exp_lower.fit,output/exp_upper.fit")'
  ```

## Reference

- Cyburt et al., ApJS 189 (2010) 240

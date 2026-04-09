./reaclib_fit -c exp.conf --input ExpRate_central.dat --output exp_central.fit --no-gui 
./reaclib_fit -c exp.conf --input ExpRate_lower.dat   --output exp_lower.fit   --no-gui 
./reaclib_fit -c exp.conf --input ExpRate_upper.dat   --output exp_upper.fit   --no-gui 

root 'macros/plot_fit_collection.C("input/ExpRate_central.dat,input/ExpRate_lower.dat,input/ExpRate_upper.dat","output/exp_central.fit,output/exp_lower.fit,output/exp_upper.fit")'

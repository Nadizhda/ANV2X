ANV2x is an adaptive numerology technique that relies on the Jakes spectrum to model the waveform design concerning the subcarrier spacing and the guard interval.
The ANV2X_scipy.py is adopted by the ANV2X.m to set the alpha and beta values.
The python file adopts the resulting Doppler, the Delay parameters, including maximum delay and rms.
These parameters are written on a text file, that is read by the python file.
Running the matlab file, will dynamically run the python file to test and evaluate the transmitted signal at each given parameters,
 including subcarrier spacing, guard interval, doppler, delay,and the transmitted symbols.
 The values of the numerology, subcarrier spacing, guard interval, alpha nd beta are then saved to a stack along with the maximum SE.
 

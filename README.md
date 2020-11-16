# Corner_reflectors
FitTest2_30350.m is a function used in the Project30350 script. Tries what polynomial of degree n best fits satellite position.
Always uses n+1 data points (i.e. one state vector adds each time when trying new degree). 

FitTest30350.m is a function used when fitting satellite position. Found degree of polynomial that best fitted all state vectors.
It was used in Project30350 script, bus FitTest2_30350.m is now used and gives better results.

Also, matlab.mat are the zero-doppler azimuth time variables obtained from method 1 (from Jesper), and they're used in the script.

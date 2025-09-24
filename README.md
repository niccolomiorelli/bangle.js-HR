Heart rate algorithms for BangleJS
==================================

Algorithms for heart rate for the Bangle JS.

Compile:

```bash
cmake --build ./build/ --target all
```

Run:

```bash
./build/heartrate ../bangle_heartrate_raw/ ../resultshr
```

Input files have the following columns: 
 ms,ms (abs),PPG,ACCx,ACCy,ACCz,GT_bangle,GT_polar

Output files returns the following columns:
 time,Dummy,Espruino,FFT,Autocorrelation,Oxford,PanTompkins,Autocorrelation2,Algo1,FFT2,GT_Bangle,GT_polar
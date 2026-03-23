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
 ms,ms (abs),PPG,ACCx,ACCy,ACCz,GT_bangle,GT_polar(,GT_cosmed)

Output files returns the following columns:
 time,time (abs),Dummy,Espruino,FFT,Autocorrelation,Oxford,PanTompkins,Autocorrelation2,FFT2,TRUST_PPG1,TRUST_PPG2,TRUST_PPG3,TRUST_PPG4,GT_Bangle,GT_polar(,GT_cosmed)
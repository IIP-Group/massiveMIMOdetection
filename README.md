# massiveMIMOdetection
Simple massive MIMO simulator that includes several data-detectors

(c) 2020 Christoph Studer and Oscar Castañeda
e-mail: studer@ethz.ch & caoscar@ethz.ch

### Important information

If you are using the simulator (or parts of it) for a publication, please consider citing our papers:

Oscar Castañeda, Tom Goldstein, and Christoph Studer, "Data Detection in Large Multi-Antenna Wireless Systems via Approximate Semidefinite Relaxation," IEEE Transactions on Circuits and Systems I: Regular Papers, vol. 63, no. 12, pp. 2334-2346, Dec. 2016.

Charles Jeon, Oscar Castañeda, and Christoph Studer, "A 354 Mb/s 0.37 mm2 151 mW 32-User 256-QAM Near-MAP Soft-Input Soft-Output Massive MU-MIMO Data Detector in 28nm CMOS," IEEE Solid-State Circuits Letters, vol. 2, no. 9, pp. 127-130, Oct. 2019.

and clearly mention this in your paper.

### How to start a simulation:

Simply run

```sh
detection_MIMO_sim
```

which starts a simulation of a 32 BS antenna, 16 user, QPSK massive MIMO system using several data-detectors. After completing, the simulator will generate the following figure:

![](doc/default_results.png?raw=true "")

You can specify your own system and simulation parameters by passing your own "par"-structure (see the simulator for an example). The simulator also includes other data-detectors that are not set to run in the predefined parameters, but you can add them by modifying "par.detector". Please note that several methods have parameters that can be tuned to improve their performance. The simulation runs in an i.i.d. Rayleigh-fading channel, but you can also use a simple line-of-sight (LoS) channel by setting "par.los" to 1.

We highly recommend you to execute the code step-by-step (using MATLAB's debug mode) in order to get a detailed understanding of the simulator.

### Data detectors in this simulator

- SIMO lower bound
- Maximum-Likelihood (ML) using sphere decoding `[Studer & Bölcskei, 2010`]
- Maximum Ratio Combining (MRC)
- unbiased Zero Forcing (ZF)
- unbiased linear Minimum Mean Squared Error (MMSE)
- Semidefinite Relaxation (SDR)
- Triangular Approximate Semidefinite Relaxation (TASER) `[Castañeda, Goldstein & Studer, 2016`]
- Row-By-Row (RBR) `[Wai, Ma & Man-Cho So, 2011`]
- Large MIMO Approximate message passing (LAMA) `[Jeon, Ghods, Maleki & Studer, 2015`]
- ADMM-based Infinity-Norm (ADMIN) `[Shahabuddin, Juntti & Studer, 2017`]
- BOX
- Optimized Coordinate Descent (OCD) `[Wu, Dick, Cavallaro & Studer, 2016`]
- LLL-Lattice-Reduction with Decision Feedback and ZF remap (LR_LLL_DFE_rZF) `[Windpassinger & Fischer, 2013`]
- K-Best

### Notes

* To use the 'SDR' data-detector, you need to first install CVX [http://cvxr.com/cvx/]

### Version history
* Version 0.2: caoscar@ethz.ch - added new detectors, support for channel estimation, and separated functions in files
* Version 0.1: caoscar@ethz.ch - initial version for GitHub release

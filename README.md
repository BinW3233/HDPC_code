# HDPC_code

This repository contains simulation code for Wang, Audette, Schneider, and Aljadeff, "Desegregation of neuronal predictive processing," bioRxiv (2024).

The code studies predictive-processing neural networks under steady-state, pulse-input, sparse-connectivity, hierarchical, and inhibitory-neuron conditions. Most files are MATLAB scripts; the inhibitory-neuron example is a Wolfram Mathematica notebook that loads a Wolfram package from this repository.

## Requirements

- MATLAB, tested with recent MATLAB releases.
- MATLAB Statistics and Machine Learning Toolbox for `mvnrnd`.
- Wolfram Mathematica 13 or later for `Inhibitory_response.nb`.

## Files

| File | Language | Purpose | Main output |
| --- | --- | --- | --- |
| `SS_network.m` | MATLAB | Simulates steady-state responses in the base predictive-processing network. | `ss_response.mat` with `hf` and `rf`. |
| `SS_TwoPair.m` | MATLAB | Simulates four match and mismatch conditions for two stimulus pairs. These responses are used to characterize cell types for Fig. 4c. | `PairedRep_*.mat` with `mmx1`, `m1`, `mmy2`, and `m2`. |
| `Pulse_Input_Network.m` | MATLAB | Simulates time-dependent responses to pulse-like `x` and `y` inputs, plus `x`-only and `y`-only controls. | `TDcode_longtxy_*.mat` with time vectors, voltage responses, weights, and parameters. |
| `SS_hirachical_net.m` | MATLAB | Simulates steady-state responses in the hierarchical predictive-processing network. The filename keeps the original spelling used by the project. | `Hierachical_response.mat` with `hf` and `rf`. |
| `Sparse_network.m` | MATLAB | Simulates time-dependent responses in a sparsified, asymmetric predictive-processing network. | `Sparse_NN.mat` with `t1`, `h1`, `r1`, weights, and parameters. |
| `IneuronRep.m` | Wolfram Language | Defines helper functions for inhibitory-neuron response calculations. | Functions loaded by `Inhibitory_response.nb`. |
| `Inhibitory_response.nb` | Mathematica notebook | Runs inhibitory-neuron response calculations by loading `IneuronRep.m`. | Notebook variables `hImis`, `hImisy`, and `hIfull`. |

## Basic Usage

Run each MATLAB script from the repository root so output files are written next to the code:

```matlab
SS_network
SS_TwoPair
Pulse_Input_Network
Sparse_network
SS_hirachical_net
```

To run the inhibitory-neuron notebook, open `Inhibitory_response.nb` in Mathematica from the repository root. The notebook uses:

```wolfram
Get["IneuronRep.m"]
```

so `IneuronRep.m` must remain in the same directory.

## Changing Stimulus Conditions

Stimulus conditions are controlled by the vectors `x` and `y` in the MATLAB scripts. They should stay as `p`-by-1 column vectors because the model multiplies them by the sampled weight matrices `w` and `v`.

Common parameters near the top of each script include:

- `nsize` or `size`: number of neurons.
- `p`: number of stimulus pairs.
- `mu`: learning or stimulus-correlation parameter.
- `b`: model gain parameter.
- `theta`: firing threshold.
- `tf`, `tpre`, `txy`, and `delta_t`: simulation timing parameters.

## Notes

- Scripts use random weight vectors, so repeated runs produce different `.mat` files unless you set a random seed with `rng(...)` before sampling.
- The scripts are written as direct, editable simulation files. Change parameters near the top of a file, run the script, and inspect the saved `.mat` output.
- Generated `.mat` files are analysis outputs and are not required to run the source code.

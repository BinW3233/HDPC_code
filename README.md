# HDPC_code
This repository contains all the code files for the preprint, Wang, Audette, Schneider, & Aljadeff. "Desegregation of neuronal predictive processing." bioRxiv (2024).


1. "SS_network.m" will generate the steady state neural responses in the predictive processing network.
2. "SS_TwoPair.m" will generate the steady state neural responses for match and mismatch conditions for two stimulus-pairs. These responses are used to generate Fig.4c in the main text and allows to characterize different function cell types. 
3. "Pulse_Input_Network.m" will generate the time-dependent neural responses for a pair of pulse-like stimuli.
4. "SS_hierachical_net.m" will generate steady state neural repsonses in the hierachical predictive processing network.
5. "Inhibitory_response.nb" will generate steady state responses for inhibitory neurons. Please make sure that "IneuronRep.m" is included in the same directory.
6. "Sparse_network.m" will generate the time-dependent neural responses for sparsified, asymmetric predictive processsing network.

Stimulus condition can be specified by changing the values of variables, x and y in the code. Other relevant parameters are also defined and modifiable in the coresponding code file. See the comments in the code file for details.

Note that "Inhibitory_response.nb" is a Mathematica code file. All other code files are Matlab codes.

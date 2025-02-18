# HDPC_code
This repository contains all the code files for the preprint, Wang, Audette, Schneider, & Aljadeff. "Desegregation of neuronal predictive processing." bioRxiv (2024).


1. "SS_network.m" will generate the steady state neural responses in the predictive processing network. 
2. "Pulse_Input_Network.m" will generate the time-dependent neural responses for a pair of pulse-like stimuli.
3. "SS_hierachical_net.m" will generate steady state neural repsonses in the hierachical predictive processing network.
4. "Inhibitory_response.nb" will generate steady state responses for inhibitory neurons. Please make sure that "IneuronRep.m" is included in the same directory.
5. "Sparse_network.m" will generate the time-dependent neural responses for sparsified, asymmetric predictive processsing network.

Stimulus condition can be specified by changing the values of variables, x and y in the code. Other relevant parameters are also defined and modifiable in the coresponding code file. See the comments in the code file for details.

Note that "Inhibitory_response.nb" is a Mathematica code file. All other code files are Matlab codes.

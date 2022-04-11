README_Aud_glm_decode.txt
A.D. Ramirez
May 11, 2012 


code for MAP decoding of spectrograms using a GLM model for the likelihood   of spiking and spectro-temporal correlations estimated from a set of stimuli


Note that this demo in both content and pedagogical style borrow from 
Johnathan Pillow's GLM point process code available at
http://pillowlab.cps.utexas.edu/code.html

INSTALLATION
============
(1) Unnpack the archive 
(2) Launch matlab and cd into the main directory containing the code
    (e.g. 'cd /GLM_Aud_Decode/').
(3) type addpath(genpath(pwd))

USE 
===
   Examine scripts in sub-directory "scripts/" to see simple
    examples of simulation/fitting/reconstruction  
Test Scripts:
-------------
1. decode_demo_WNprior.m - reconstructs correlated spectrogram using a 
			   standard normal prior and simulated GLM  
                           responses

2. decode_demo_STprior.m - reconstructs correlated spectrogram using a 
			   correlated normal prior and simulated GLM  
                           responses

3. decodeHandel_demo_STprior.m - 
reconstructs (correlated) spectrogram corresponding to Handel's Hallelujah chorus using a correlated normal prior and simulated GLM  responses


NOTES
=====
- in the first two demonstrations time is represented in units of "stimulus frames", which is controlled by the global variable RefreshRate (Hz).  Thus, for
example, if RefreshRate=100, each unit of time is 10ms, and a
post-spike kernel discretized in time bins of width dt=.02 has time
bins of length 0.2 ms.

- fitting code relies on the matlab optimization toolbox ("fminunc",
"fmincon").

- this code is published under the GNU General Public License.  You
may freely redistribute it and/or modify it under the terms of the GNU
General Public License (GPLv3), available at
http://www.opensource.org/licenses/gpl-3.0.html
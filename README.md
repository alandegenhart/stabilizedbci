# Stabilized BCI

## Welcome!

This project contains companion code for BCI stabilization and example data for the paper "A stabilized brain-computer interface based on neural manifold alignment" by A. D. Degenhart, W. E. Bishop, E. R. Oby, E. C. Tyler-Kabara, S. M. Chase, A. P. Batista and B. M. Yu published in Nature Biomedical Engineering.

## Getting started with stabilization

The code for this project has been organized in a modular way to make it straightforward to apply stabilization to existing BCI decoders and to also demonstrate how stabilization was used with a Kalman filter, as reported in the original paper. The best place to get started is with the script 'stabilizationExample.m', which provides a 'bare-bones' example of how to apply stabilization to neural data.  This script is intended to help researchers who have an existing decoder they simply want to apply stabilization to.  The script 'stabilizedKalmanFilterExample.m' is provided to reproduce the closed-loop online results that were achieved with a stabilized Kalman filter, as reported in the original paper.  All the code to train a base Kalman filter and apply stabilization is demonstrated in this script.

## Organization of the code
In addition to the two top-level example scripts, the core code is contained in 'code' folder, and is organized as follow:

- stabilization: This subfolder contains the core functions needed to apply stabilization to your own decoder.  The three main functions needed are 'fitBaseStabilizer.m' for fitting an initial stabilizer to neural data collected during a calibration session, 'updateStabilizer.m' for updating the stabilizer with neural data collected during normal BCI use and 'getStabilizationMatrices.m' for getting a set of matrices which can be applied to linearly extract stabilized neural signals.  The use of these functions are all demonstrated in the 'stabilizationExample.m' script mentioned above. 

- stabilizedKalman: This folder contains code specific to using a stabilized Kalman filter.  It builds on the code in the 'stabilization' folder.  The two main functions are 'stabilizedKalmanFit.m' for fitting a baseline Kalman filter and stabilizer to neural data collected during a calibration session and 'stabilizedKalmanDecodeMatrixForm.m' for decoding with a stabilized Kalman filter.  For applying stabilization updates, the 'updateStabilizer.m' function (under the stabilization folder) can be directly applied.  The use of these functions are demonstrated in the 'stabilizedKalmanFilterExample.m' mentioned above.

- utilities: These are basic support functions. 

## Description of the data

Data from an example experiment (performed on 3/25/206) is saved in the data folder as a .mat file.  

Note that in the code and documentation below an instability is often referred to as a "perturbation," and for this work the two terms can be used interchangeably. 

The .mat file containing the data contains the following fields: 

- trialData: a structure of length T, where T is the number of trials.  It will have the following fields: 

	- id: The id of each trial.

	- success: True if the trial was a success.

	- binCounts: The binned counts for the trial of shape electrodes*time_steps; these are the raw counts as recorded from all electrodes on the array; note that not all of these electrodes were good, so some of them were not used (see goodElectrodes below). 

	- refPos: The position of the cursor at the start of the bin.

	- tgtPos: The position of the target.

	- decodedVl: The decoded velocity.

	- moveBin: The index of bin in which the cursor was first allowed to move.

	- decoderInd: The index of the original saved decoder used to decode the trial in the closed-loop work.  Each stabilized-decoder was recorded as a new decoder, with it's own index. 
     
	- trainKinematics: kinematics used for initial calibration.

	- type: The type of trial, which will be one of the following strings:
        	
		- observation: Initial trials in which the monkey watched but had no control over the movement of the cursor; these are the first calibration trials. 
		- gradualTrian: Trials when we were initially training the decoder in a supervised manner; these are still calibration trials. 
		- baselineEvaluation: Trials when the fully trained intuitive decoder was used and no perturbation was applied. This is the start of normal BCI use. 
		- perturbationEvaluation: trials where the baseline decoder was used to decode and a perturbation was applied.
		- stabilizationPerturbation: the first trials after we have applied the perturbation and stabilization is running but the stabilizer has not been updated yet.
		- stabilization: trials that were decoded with an updated stabilizer (after the perturbation has been applied).
		- stabilizationEvaluation: trials using the final stabilizer in the face of a perturbation.

- goodElectrodes: Indices into trialData.binCounts of good electrodes

- pertTrials: An array of indices into trialData, indicating which trials the perturbation was applied to.  

- pertParams: A structure with the parameters for the perturbation. 

- baselineTrainParams - Settings used for training the baseline decoder.

- stabilizationParams - Settings used for performing stabilization updates.

- decoderTrainInfo: A structure with information that is needed if we are to re-train the baseline decoder or re-run stabilizer to get the same decoder/stabilizers that were used in the online, closed loop experiments. Each stabilization update was considered, for the purposes of the MATLAB code, a new decoder. decododerTrainInfo(i) contains the essential information necessary to reproduce each decoder/stabilizer.  decoderTrainInfo(1) contains information for the baseline Kalman filter and stabilizer. decoderTrainInfo(2...) contains information for the stabilizer updates.  This structure will have the following fields:
	- trainTrials: a structure with a length equal to the number of trials used for training/stabilization.  It will contain the fields 'trialNum' with the index of a trial used for training/stabilization and 'trainBins' with indices of bins used for training/stabilization from that trial. 
	- cInit: The initial loading matrix used to initialize EM for fitting the stabilizer
	- psiInit: The initial matrix of private variances used to initialize EM for fitting the stabilizer

## Getting help

For questions, please contact William Bishop at <wbishop@janelia.hhmi.org> or Alan Degenhart at <add19@pitt.edu>.


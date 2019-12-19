% This is a "bare bones" example of how to perform neural stabilization.
% Stabilization can be combined with many different decoders, so this code
% only shows how to:
%
%   1) Fit an initial stabilizer to neural data (this would typically be 
%   done using neural data collected during initial BCI calibration)
%
%   2) Update the stabilizer with newly collected neural data (this would
%   typically be done using neural data collected during normal BCI use)
%
%   3) Get stabilized signals out of the stabilizer
%
% Once stabilized signals are extracted from the stabilizer, these can be
% fed into a neural decoder.  To see an example of how stabilization can 
% be combined with a Kalman filter, see stabilizedKalmanFilterExample.m
%
% Note: The code below assumes you are in the main project directory (the
% directory containing this script) when running it. 

%%  =======================================================================
% Add the required paths to use the stabilization code
startDir = pwd;
cd('code'); 
addStabilizationProjectPaths();
cd(startDir); 

%%  =======================================================================
%   We load some neural data to work with.  We will work with neural
%   data that is nominally stable (i.e., collected in a single day while
%   the monkey was working a baseline decoder) and we will then add some
%   instabilities to it to demonstrate how stabilization works.  
%
%   When doing this we also downselect to only the good electrodes on the
%   array (3 were bad). 
%
%   The result is that each entry of neuralData will contain nominally 
%   stable counts for 93 electrodes and for one trial.  

dataFile = fullfile('data', '20160325');
data = load(dataFile); 
baselineTrials = find(strcmp({data.trialData.type}, 'baselineEvaluation'));

goodElectrodes = data.goodElectrodes;
nGoodElectrodes = length(goodElectrodes);
neuralData = cellfun(@(x) x(:, goodElectrodes)', ...
    {data.trialData(baselineTrials).binCounts}, 'UniformOutput', false);

%%  =======================================================================
%   We now split up our neural data into calibration trials (trials we fit 
%   the initial stabilizer to) and instability trials (trials where we will
%   apply a simulated instability).

nCalibrationTrials = 100; % Can adjust this to change the number of calibration trials
calibrationData = neuralData(1:nCalibrationTrials); 
instabilityData = neuralData(nCalibrationTrials+1:end); 

%%  =======================================================================
%   We now create the instability. What we will do is simulate starting
%   with a set of 83 electrodes for calibration.  We will then apply a 
%   simulated instability to these electrodes.  In particular we will:
%
%       1) Change the tuning of the first 10 electrodes (by swapping their
%       activity with the other set of 10 electrodes that were not included
%       in the original 83). 
%
%       2) Silencing the next 5 electrodes
%
%       3) Adding an offset in spike counts to the remaining electrodes

% These parameters can be adjusted to change the type of instability
nSwapElectrodes = 10; 
nSilenceElectrodes = 5;
nOffsetElectrodes = 68; 

% State how many electrodes we want to simulate recording from 
nBaseElectrodes = 83;

% Make sure we have enough electrodes to apply the requested instability
if  2*nSwapElectrodes + nSilenceElectrodes + nOffsetElectrodes > nGoodElectrodes
    error('Not enough original electrodes to apply the requested instability.')
end

% Make sure we have enough base electrodes to apply the requested instability
if  nSwapElectrodes + nSilenceElectrodes + nOffsetElectrodes > nBaseElectrodes
    error('Number of base electrodes needs to be increased to accomodate the requested instability.')
end

% Before applying instabilities, we make a copy of the original counts for 
% the instability trials; this will allow us to form a "ground truth" of 
% what should be coming out of the stabilizer if we have no instabilies 
% and don't update the stabilizer. 
origInstabilityData = cellfun(@(x) x(1:nBaseElectrodes,:), instabilityData, 'UniformOutput', false); 

% Calcualte the offsets we will apply - we introduce a bias here to make
% the effect of the instability more pronouced (as random offsets across
% neurons could counteract one another).
randomOffsets = randn(nOffsetElectrodes, 1) + 2; 

% Apply the offsets here
for tI = 1:length(instabilityData)
    % Tuning changes
    instabilityData{tI}(1:nSwapElectrodes,:) = instabilityData{tI}(nGoodElectrodes-nSwapElectrodes+1:end,:);  
    % Silence electrodes
    instabilityData{tI}(nSwapElectrodes+1:nSwapElectrodes+nSilenceElectrodes,:) = 0;
    % Apply random offsets
    instabilityData{tI}(nSwapElectrodes+nSilenceElectrodes+1:nSwapElectrodes+nSilenceElectrodes+nOffsetElectrodes,:) = ...
        bsxfun(@plus, instabilityData{tI}(nSwapElectrodes+nSilenceElectrodes+1:nSwapElectrodes+nSilenceElectrodes+nOffsetElectrodes,:), randomOffsets); 
end

% Down-select to the base electrodes
calibrationData = cellfun(@(x) x(1:nBaseElectrodes, :), calibrationData, 'UniformOutput', false); 
instabilityData = cellfun(@(x) x(1:nBaseElectrodes, :), instabilityData, 'UniformOutput', false); 

%%  =======================================================================
%   We now fit the baseline stabilizer - this is the stabilizer that
%   would be fit during normal BCI calibration. 

nLatents = 10; % Number of latents to use for stabilizaton; see methods and
               % supplemental figure 7 of the paper for a discusson on 
               % choosing this value

baselineStabilizer = fitBaseStabilizer(calibrationData, nLatents);

%%  =======================================================================
%   We now update the stabilizer using a fixed number of trials during 
%   the instability period.  In practice, a buffer of neural data would
%   be maintained and the stabilizer would be updated periodically with
%   neural data in this buffer

nInstabilityUpdateTrials = 100; % This can be adjusted
stabilizerUpdateData = instabilityData(1:nInstabilityUpdateTrials); 

% Here we update the stabilizer; the important parameter is ALIGN_N which
% is the number of stable electrodes we expect in the neural data.  Note
% that for the purposes of alignment, electrodes which have only suffered
% a baseline shift, can still be considered "stable" for the purposes of
% alignment. See methods and supplemental figure 7 of the original paper 
% for guidance on choosing this number.  

% The value of ALIGN_TH less important.  This is a threshold on the l_2
% norm of an electrode's row in a loading matrix.  Electrodes with rows
% in either the loading matrix for the original or updated stabilizer with
% norms below this are immediately considered 'unstable.' We suggest
% setting this to a small value below the l_2 norm seen for normal,
% non-silent electrodes without clear artifacts in your data. 

updatedStabilizer = updateStabilizer(baselineStabilizer, stabilizerUpdateData, ...
                                     'ALIGN_N', 60, 'ALIGN_TH', .01); 

%%  =======================================================================
% We also pull out trials during the instability period where we will 
% evaluate the performance of the updated stabilizer; these are trials
% that were not used for calibration nor were they used for updaing the
% stabilizer.  We refer to these as "evaluation" trials. 

instabilityEvalData = instabilityData(nInstabilityUpdateTrials:end); 
origInstabilityEvalData = origInstabilityData(nInstabilityUpdateTrials:end); 

%%  =======================================================================
%   We now extract stabilized signals in the evaluation trials during the
%   instability period.  The output of the stabilizer is referred to as 
%   'latent state."  We will extract latent state in three ways:
%
%   1) We will extract latent state from the neural data with instabilities
%      applied using the original baseline stabilizer
%
%   2) We will extract latent state from the neural data with the 
%      instabilities applied with the updated stabilizer
%
%   3) Will extract latent state from the neural data without instabilities
%      using the original baseline stabilizer; this is to give us a
%      reference to "best possibe" performance - that which we would get
%      if there were no instabilities applied and we just kept using the
%      baseline stabilizer

% We first get the stabilizers in matrix form; this form is useful for
% online decoding
[baselineBeta, baselineO] = getStabilizatonMatrices(baselineStabilizer);
[updatedBeta, updatedO] = getStabilizatonMatrices(updatedStabilizer); 

% Now we extract latent state in all 3 ways
nEvalTrials = length(instabilityEvalData); 
baselineInstabilitiesLatentState = cell(1, nEvalTrials); 
updatedInstabilitiesLatentState = cell(1, nEvalTrials); 
baselineOriginalLatentState = cell(1, nEvalTrials); 
for tI = 1:nEvalTrials
    baselineInstabilitiesLatentState{tI} = bsxfun(@plus, baselineBeta*instabilityEvalData{tI}, baselineO);
    updatedInstabilitiesLatentState{tI} = bsxfun(@plus, updatedBeta*instabilityEvalData{tI}, updatedO);
    baselineOriginalLatentState{tI} = bsxfun(@plus, baselineBeta*origInstabilityEvalData{tI}, baselineO);
end

%%  =======================================================================
%   Now we look at an example trial; we will look at what the latent
%   state would have been had we continued to use the baseline stabilizer
%   without instabilities and we will see how close we come to recovering 
%   this when apply instabilities but apply stabilization.  We also see
%   what happens had we continued to use the baseline stabilizer without
%   updating it. 

exTrial = 10;
figure();
for lI = 1:nLatents
    subplot(ceil(nLatents/2)+1, 2, lI); 
    plot(baselineOriginalLatentState{exTrial}(lI, :), 'ko-'); 
    hold on;
    plot(updatedInstabilitiesLatentState{exTrial}(lI, :), 'b.-'); 
    plot(baselineInstabilitiesLatentState{exTrial}(lI, :), 'r--');
    ylabel(['Latent ', num2str(lI), ' (a.u.)']); 
    if lI == nLatents
        xlabel('Bin Number')
        leg = legend('Baseline stabilizer with original data', 'Updated stabilizer with instabiliites', 'Baseline stabilizer with instabilities');
        set(leg, 'Position', [.65, .05, .2, .1], 'FontSize', 15); 
    end
    set(gcf, 'Position', [0 0, 1000, 800]); 
end














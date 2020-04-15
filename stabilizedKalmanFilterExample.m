% This script demonstrates how to use the code provided in this project to
% stabilize a Kalman Filter. 
%
% This script will reproduce results of a single-day closed-loop
% experiment, as described in the original paper.  
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
%   We load the neural data as well as parameters we need for reproducing
%   the closed-loop results here.  

dataFile = fullfile('data', '20160325');
data = load(dataFile);
trialData = data.trialData;
pertTrials = data.pertTrials;
pertParams = data.pertParams;
decoderTrainInfo = data.decoderTrainInfo;
baselineTrainParams = data.baselineTrainParams;
stabilizationParams = data.stabilizationParams;

%% ========================================================================
%  Apply perturbation to the counts; we apply the same simulated
%  perturbations that were applied in the closed-loop experiment
%  ========================================================================
trialData = perturbCounts(trialData, pertTrials, pertParams);

%% ========================================================================
%  Train the base decoder; in the closed-loop experiment an iterative
%  calibration procedure was used, in which multiple decoders were fit.
%  Here we just fit the final baseline decoder that was used in the rest
%  of the experiment.  
%  ========================================================================

N_LATENTS = 10; % A dimensionality of 10 was used for the manifold

% Find the trials used to train the baseline decoder
baselineDecoderInd = trialData(find(strcmp({trialData.type}, 'baselineEvaluation'), 1)).decoderInd;
baselineDecoderTrainTrials= decoderTrainInfo(baselineDecoderInd).trainTrials;

% Get the data used to train the baseline decoder
nTrainTrials = length(baselineDecoderTrainTrials);
x = cell(1, nTrainTrials); % Will hold training kinematic data
y = cell(1, nTrainTrials); % Will hold training neural data
for tI = 1:nTrainTrials
    curTrialNum = baselineDecoderTrainTrials(tI).trialNum;
    curBins = baselineDecoderTrainTrials(tI).trainBins;
    x{tI} = trialData(curTrialNum).trainKinematics(curBins, :)';
    y{tI} = trialData(curTrialNum).pertBinCounts(curBins,pertParams.baseElectrodes)';
end

% Train the baseline decoder - note in practice C_INIT and PSI_INIT would
% not be used; here we use them to ensure we start the EM algorithm for
% fitting stabilizer from the same initial conditions they were started
% from in the closed-loop work; this is to ensure reproducability in 
% this script, but in practice these optional inputs can be omitted,
% allowing EM to start from random initial conditions
%
baselineDecoder = stabilizedKalmanFit(x, y, N_LATENTS, ...
    'N_FA_RESTARTS', baselineTrainParams.N_FA_RESTARTS, ...
    'LL_DIFF_THRESH', baselineTrainParams.LL_DIFF_THRESH, ...
    'MAX_N_ITS', baselineTrainParams.MAX_N_ITS, ...
    'MIN_PRIV_VAR', baselineTrainParams.MIN_PRIV_VAR, ...
    'FIT_STATE_NOISE', baselineTrainParams.FIT_STATE_NOISE, ...
    'STATE_NOISE_VAR', baselineTrainParams.STATE_NOISE_VAR, ...
    'C_INIT', decoderTrainInfo(baselineDecoderInd).cInit, ...
    'PSI_INIT', decoderTrainInfo(baselineDecoderInd).psiInit);

%% ========================================================================
%  Run all stabilization updates - in practice, this would be done online,
%  running stabilization updates as more data is collected
%  ========================================================================

% Find indices of all stabilizer updates - these indices are needed to
% allow us to determine which trials were used for each stabilizer update
stabilizationUpdateTrials = [find(strcmp({trialData.type}, 'stabilization')), ...
                             find(strcmp({trialData.type}, 'stabilizationEvaluation'))];
stabilizationUpdateInds = unique([trialData(stabilizationUpdateTrials).decoderInd]); 
nStabilizerUpdates = length(stabilizationUpdateInds);
stabilizedDecoders = cell(1, nStabilizerUpdates);

% Perform the stabilizer updates
for sI = 1:nStabilizerUpdates
    
    % Pull out the neural data that was used for this stabilizer update
    stabilizerTrialInfo = decoderTrainInfo(stabilizationUpdateInds(sI)).trainTrials;
    nStabilizationTrials = length(stabilizerTrialInfo);
    y = cell(1, nStabilizationTrials);
    for tI = 1:nStabilizationTrials
        curTrialNum = stabilizerTrialInfo(tI).trialNum;
        curBins = stabilizerTrialInfo(tI).trainBins;
        y{tI} = trialData(curTrialNum).pertBinCounts(curBins,pertParams.baseElectrodes)';
    end
    
    % Perform stabilization; just as with fitting the baseline decoder
    % and stabilizer, in practice, the optional inputs C_INIT and PSI_INIT
    % can be omitted.
    stabilizedDecoders{sI}.kf = baselineDecoder.kf;
    stabilizedDecoders{sI}.fa = updateStabilizer(baselineDecoder.fa,  y, ...
        'N_FA_RESTARTS', stabilizationParams.N_FA_RESTARTS, ...
        'LL_DIFF_THRESH', stabilizationParams.LL_DIFF_THRESH, ...
        'MAX_N_ITS', stabilizationParams.MAX_N_ITS, ...
        'ALIGN_TH', stabilizationParams.ALIGN_TH, ...
        'ALIGN_N', stabilizationParams.ALIGN_N, ...
        'MIN_PRIV_VAR', stabilizationParams.MIN_PRIV_VAR, ...
        'C_INIT', decoderTrainInfo(stabilizationUpdateInds(sI)).cInit, ...
        'PSI_INIT', decoderTrainInfo(stabilizationUpdateInds(sI)).psiInit);
end

%% ========================================================================
%  Put the baseline and each stabilized decoder into matrix form - this is
%  possible using the steady state representation of the Kalman Filter;
%
%  NOTE: As described in the methods of the paper, the parameters of the
%  Kalman filter were calculated in a manner slightly different from the 
%  standard equations for the steady-state Kalman filter.  To calculate these
%  parameters in the same manner as they were in the closed loop work, we
%  have set 'COMPUTE_METHOD' below to 'OLD'.  However, if using this
%  code in future work, 'COMPUTE_METHOD' should be set to 'CORRECT' :)
%  ========================================================================
baselineDecoderMatForm = convertStabilizedKalmanToMatrixForm(baselineDecoder, 'COMPUTE_METHOD', 'OLD');
stabilizedDecodersMatForm = cell(1, nStabilizerUpdates); 
for sI = 1:nStabilizerUpdates
    stabilizedDecodersMatForm{sI} = convertStabilizedKalmanToMatrixForm(stabilizedDecoders{sI}, 'COMPUTE_METHOD', 'OLD');
end

%% ========================================================================
% Now we re-decode all trials from the original experiment in which
% the baseline decoder was used with and without stabilization
%  ========================================================================

% Find the index indicating when the baseline decoder was applied with
% perturbed neural data
baselineDecoderWithPertInd = trialData(find(strcmp({trialData.type}, 'perturbationEvaluation'), 1)).decoderInd;

% Note the decoder indices when the baseline decoder was applied w/ and w/o
% stabilization
baselineDecoderInds = [baselineDecoderInd, baselineDecoderWithPertInd];

% Note all decoder indices for trials we will decode again - we will look
% for trials which were originally decoded with a decoder with one of 
% these indices and re-decode them; these correspond to all trials which 
% were decoded with the baseline decoder w/ and w/o stabilization 
redecodeInds = [baselineDecoderInds, stabilizationUpdateInds];

nOrigTrials = length(trialData); % Number of original trials, including calibration trials 
decodedVelocities = cell(1, nOrigTrials); 

for tI = 1:nOrigTrials
    if any(trialData(tI).decoderInd == redecodeInds)
        
        if any(trialData(tI).decoderInd == baselineDecoderInds)
            decodedVelocities(tI) = stabilizedKalmanDecodeMatrixForm(baselineDecoderMatForm, ...
                {trialData(tI).pertBinCounts(:, pertParams.baseElectrodes)'});
        else
            stabilizedInd = find(stabilizationUpdateInds == trialData(tI).decoderInd);
            decodedVelocities(tI) = stabilizedKalmanDecodeMatrixForm(stabilizedDecodersMatForm{stabilizedInd}, ...
                {trialData(tI).pertBinCounts(:, pertParams.baseElectrodes)'});
        end
    end
end

%% ========================================================================
% Now we go through and make sure the decode is close to the original for
% all decoded trials
%  ========================================================================
trialDiffs = cell(1, nOrigTrials); 
for tI = 1:nOrigTrials
    if any(trialData(tI).decoderInd == redecodeInds)
        trialDiffs{tI} = trialData(tI).decodedVl' - decodedVelocities{tI};
    end
end
trialErs = cellfun(@(x) mean(mean(x.^2)), trialDiffs);
figure;
plot(trialErs)
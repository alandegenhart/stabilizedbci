% processTrialData  Process trial data into blocks for subsequent analysis.
%
% This function formats the provided session trial data into blocks
% corresponding to distinct phases of an experimental session.
%
% Usage:
%   [D] = processTrialData(trialData)
%
% Inputs:
%   trialData   Structure array of data for individual trials
%
% Outputs:
%   D           Structure array of data segmented into blocks
%
% @ Alan Degenhart -- alan.degenhart@gmail.com

function [D] = processTrialData(trialData)

% Constant parameters
blockSize = 16;
baselineDecoderNum = 10;

% Initialize matrices
totalTrials = length(trialData);
successCode = [trialData.success];
onsetTime = [trialData.cursorOnsetTime];
offsetTime = [trialData.trialEndTime];
acquireTime = double(offsetTime - onsetTime);

% Identify decoder blocks. The decoder used was changed over the course of
% the experiment. First identify contiguous blocks of trials with the same
% decoder. Note that the decoder ID does not always increase -- for
% baseline evaluation/washout blocks at the end of a session the decoder
% was switched back to the baseline decoder (without any applied
% instabilities).
decoderNum = [trialData.decoderInd];
decoderNum(isnan(decoderNum)) = -1;  % Set ID for trials w/o a decoder

% Identify post-baseline evaluation blocks. Trials were briefly stopped
% after the baseline evaluation block for some sessions in order to select
% an instability. Following selection of the instabiliy, ~32 trials were
% run with the baseline decoder so that the animal was engaged in the task
% when the instability was first applied. To allow these trials to be
% separated out, give them a different decoder number.
preInstDecoderNum = 100;
baselineTrials = diff(decoderNum == baselineDecoderNum);
baselineOnset = find(baselineTrials == 1, 1, 'first') + 1;
baselineOffset = find(baselineTrials == -1, 1, 'first');
decoderNum(baselineOnset + 128:baselineOffset) = preInstDecoderNum;

% Find changes in decoder number. These indicate different 'blocks' of
% trials to segment. Also reset the decoder number for the 'pre-instability
% baseline' trials, as the actual decoder number needs to be correct for
% subsequent analyses.
blockOnsetIdx = find(diff(decoderNum) ~= 0) + 1;
decoderNum(decoderNum == preInstDecoderNum) = baselineDecoderNum;
blockDecoderInd = decoderNum(blockOnsetIdx);

% Remove any blocks during decoder calibration. Decoder IDs 1-9 correspond
% to the intermediate decoders created during the "gradual training" block
% of trials, and should not be analyzed.
mask = blockDecoderInd >= 10;
blockOnsetIdx = blockOnsetIdx(mask);
blockDecoderInd = blockDecoderInd(mask);
nDecoderBlocks = length(blockDecoderInd);

% Iterate over blocks and get data for each trial. 
D = repmat(defineBlockStructure(), nDecoderBlocks, 1);
TDcell = cell(nDecoderBlocks,1);
for i = 1:nDecoderBlocks
    % Get onset and offset for decoder
    onset = blockOnsetIdx(i);
    if i == nDecoderBlocks
        offset = totalTrials;
    else
        offset = blockOnsetIdx(i+1) - 1;
    end
    
    % Get success code and acquisition time for current decoder
    decSuccessCode = successCode(onset:offset);
    decAcquireTime = acquireTime(onset:offset);
    decTrials = length(decSuccessCode);
    
    % Loop over blocks and calculate success rate.  Note: if there is at
    % least one full block, any remaining trials are discarded.  Thus, the
    % 'block' metrics should not be used when analyzing data from the
    % stitching block, where the decoder number changes every ~16 trials.
    % Instead, the blockwise metrics should be calculated for all of the
    % trials with a specific decoder number.
    blockOnset = 1:blockSize:decTrials;
    blockOffset = blockOnset + blockSize - 1;
    if blockOffset(end) > decTrials
        blockOnset = blockOnset(1:(end-1));
        blockOffset = blockOffset(1:(end-1));
    end
    nBlocks = length(blockOnset);
    
    % If there is not at least one block, still compute sucess rate over
    % all available trials.
    if nBlocks == 0
        nBlocks = 1;
        blockOnset = 1;
        blockOffset = decTrials;
        warning('Fewer than %d trials exist for decoder %d.', ...
            blockSize,decoderNum(onset))
    end
    
    blockSuccessRate = nan(1,nBlocks);
    blockAcquireTimeAllMean = nan(1,nBlocks);
    blockAcquireTimeSuccessMean = nan(1,nBlocks);
    blockAcquireTimeMedian = nan(1,nBlocks);
    blockRewardRate = nan(1,nBlocks);
    for j = 1:nBlocks
        % Get success rate and acquisition time data for current block
        blockInds = blockOnset(j):blockOffset(j);
        tempSuccessRate = decSuccessCode(blockInds);
        tempAcquireTimeAll = decAcquireTime(blockInds);
        tempAcquireTimeSuccess = tempAcquireTimeAll(logical(tempSuccessRate));
        
        % Calculate success rate and mean/median acquisition time.
        % Acquisition time statistics are only calculated over successful
        % trials.
        blockSuccessRate(j) = mean(tempSuccessRate);
        blockRewardRate(j) = sum(tempSuccessRate)/sum(tempAcquireTimeAll);

        blockAcquireTimeAllMean(j) = mean(tempAcquireTimeAll);
        blockAcquireTimeSuccessMean(j) = mean(tempAcquireTimeSuccess);
        blockAcquireTimeMedian(j) = median(tempAcquireTimeAll);
    end
    
    % Put trajectory data objects into a cell
    TDcell{i} = trialData(onset:offset);
    
    % Put decoder info into structure array
    D(i).decoderNum = decoderNum(onset);
    D(i).nTrials = length(decSuccessCode);
    D(i).TD = trialData(onset:offset);
    D(i).successCode = decSuccessCode;
    D(i).acquireTime = decAcquireTime;
    D(i).successRate = sum(decSuccessCode)/D(i).nTrials;
    D(i).acquireTime_mean_all = mean(decAcquireTime);
    D(i).acquireTime_mean_suc = mean(decAcquireTime(logical(decSuccessCode)));
    D(i).rewardRate = sum(decSuccessCode)/sum(decAcquireTime); % Rewards/ms
    D(i).nBlocks = nBlocks;
    D(i).blockSize = blockSize;
    D(i).blockSuccessRate = blockSuccessRate;
    D(i).blockAcquireTimeAllMean = blockAcquireTimeAllMean;
    D(i).blockAcquireTimeSuccessMean = blockAcquireTimeSuccessMean;
    D(i).blockAcquireTimeMedian = blockAcquireTimeMedian;
    D(i).blockRewardRate = blockRewardRate;
end

% Find evaluation blocks
decoderNum = [D.decoderNum];

% Define key decoder numbers and the corresponding block indices. 
perturbedDecoderNum = baselineDecoderNum + 1;
stitchedDecoderNum = max(decoderNum);

% Find block indicies.  This is determined by the order of the blocks in
% the experiment structure. In the case of the representative experiment,
% the order was: (1) baseline, (2) pre-instability baseline, 
% (3) stabilization, (4) stabilization evaluation, 
% (5) instability evaluation, (6) baseline washout
baselineEvalInd = find(decoderNum == baselineDecoderNum, 1, 'first');
stitchingEvalInd = find(decoderNum == stitchedDecoderNum, 1, 'last');
perturbationEvalInd = find(decoderNum == perturbedDecoderNum, 1, 'last');
baselineWashoutInd = find(decoderNum == baselineDecoderNum, 1, 'last');
preInstabilityInd = find(decoderNum == preInstDecoderNum);
stitchOnset = find(decoderNum == perturbedDecoderNum, 1, 'first');
stitchOffset = stitchingEvalInd - 1;

% Set notes, evaluation, and highlight flags for baseline eval block
if ~isempty(baselineEvalInd)
    D(baselineEvalInd).blockNotes = 'Baseline Evaluation';
    D(baselineEvalInd).evaluationBlock = true;
    D(baselineEvalInd).highlightBlock = true;
end

% Set notes, evaluation, and highlight flags for stabilization eval block
if ~isempty(stitchingEvalInd)
    D(stitchingEvalInd).blockNotes = 'Stitching Evaluation';
    D(stitchingEvalInd).evaluationBlock = true;
    D(stitchingEvalInd).highlightBlock = true;
end

% Set notes, evaluation, and highlight flags for instability eval block
if ~isempty(perturbationEvalInd)
    D(perturbationEvalInd).blockNotes = 'Perturbation Evaluation';
    D(perturbationEvalInd).evaluationBlock = true;
    D(perturbationEvalInd).highlightBlock = true;
end

% Set notes, evaluation, and highlight flags for baseline washout block
if ~isempty(baselineWashoutInd)
    D(baselineWashoutInd).blockNotes = 'Baseline Washout';
    D(baselineWashoutInd).evaluationBlock = true;
    D(baselineWashoutInd).highlightBlock = true;
end

% Set notes, evaluation, and highlight flags for pre-instability block
if ~isempty(preInstabilityInd)
    D(preInstabilityInd).blockNotes = 'Pre-instability Baseline';
    D(preInstabilityInd).highlightBlock = false;
    D(preInstabilityInd).decoderNum = baselineDecoderNum;
end

% Set notes, evaluation, and highlight flags for stabilizer block
if ~isempty(stitchOffset)
    for i = stitchOnset:stitchOffset
        % Check to make sure that block information has not yet been filled
        % out.  If so, it is an indication the onset and offset were
        % incorrectly identified.
        if ~isempty(D(i).blockNotes)
            warning('Block information is not empty for decoder %d.', ...
                D(i).decoderNum)
            continue
        end
        D(i).blockNotes = 'Stitching';
        D(i).evaluationBlock = false;
        D(i).highlightBlock = false;
    end
end

end


function [D] = defineBlockStructure()

% Put decoder info into structure array
D.decoderNum = [];
D.nTrials = [];
D.TD = [];
D.targetCode = [];
D.successCode = [];
D.acquireTime = [];
D.successRate = [];
D.acquireTime_mean_all = [];
D.acquireTime_mean_suc = [];
D.rewardRate = []; % Rewards/ms
D.nBlocks = [];
D.blockSize = [];
D.blockSuccessRate = [];
D.blockAcquireTimeAllMean = [];
D.blockAcquireTimeSuccessMean = [];
D.blockAcquireTimeMedian = [];
D.blockRewardRate = [];
D.blockNotes = [];
D.evaluationBlock = false;  % Block is an evaluation block
D.highlightBlock = false;  % Use background color when plotting to highlight

end
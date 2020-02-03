% processTrialData Process trial data into blocks for subsequent analysis.
%
% 

function [E,D] = processTrialData(trialData,varargin)

% Optional arguments
saveFigs = false;
saveLoc = [];
saveData = false;
baselineDecoderNum = [];
plotPerformance = true;
plotTrajectories = true;
calcStatistics = true;
plotStatistics = true;
acqTimeOption = 'successOnly';
blockSize = 16;

% Decoder numbers
baselineDecoderNum = 10;


% Parse optional arguments
assignOpts(varargin);

% Initialize matrices
totalTrials = length(trialData);
successCode = [trialData.success];
acquireTime = ones(totalTrials,1) * 0.5;  % NOTE: Placeholder data

% Convert target positions to a unique set of target codes. These are
% useful for quickly grabbing all trials to a specific target.
% NOTE: might not need to do this if it's easier just to plot targets by
% unique positions.
targetCode = ones(totalTrials, 1);  % NOTE: Placeholder data

% Identify decoder blocks. The decoder used was changed over the course of
% the experiment. First identify contiguous blocks of trials with the same
% decoder. Note that the decoder ID does not always increase -- for
% baseline evaluation/washout blocks at the end of a session the decoder
% was switched back to the baseline decoder (without any applied
% instabilities).
decoderNum = [trialData.decoderInd];
decoderNum(isnan(decoderNum)) = -1;  % Set ID for trials w/o a decoder
blockOnsetIdx = find(diff(decoderNum) ~= 0) + 1;
blockDecoderInd = decoderNum(blockOnsetIdx);

% Remove any blocks during decoder calibration. Decoder IDs 1-9 correspond
% to the intermediate decoders created during the "gradual training" block
% of trials, and should not be analyzed.
mask = blockDecoderInd >= 10;
blockOnsetIdx = blockOnsetIdx(mask);
blockDecoderInd = blockDecoderInd(mask);
nDecoderBlocks = length(blockDecoderInd);

% NOTE: One issue with the current data format is that the data is no
% longer separated by directories. This was used for determining task
% structure somewhat. The main issue here are the last 32 trials of with
% the baseline decoder prior to the instability being introduced. In the
% current data file, there is no way to tell when the end of the baseline
% evaluation block was purely based on the trial data. In this case,
% hard-code the baseline evaluation block to be the first 128 trials with
% the baseline decoder.

% TODO: temporarily block out the pre-instability baseline trials. Can set
% these to have a different decoder ID, which should ensure that the below
% code segments them appropriately.
preInstDecoderNum = 100;

% Iterate over blocks and get data for each trial. 
D = repmat(defineBlockStructure(), nDecoderBlocks, 1);
TDcell = cell(nDecoderBlocks,1);
for i = 1:nDecoderBlocks
    % Get onset and offset for decoder
    onset = blockOnsetIdx(i);
    if i == nDecoderBlocks
        offset = totalTrials;
    else
        offset = blockOnsetIdx(i+1);
    end
    
    % Get success code and acquisition time for current decoder
    decSuccessCode = successCode(onset:offset);
    decAcquireTime = acquireTime(onset:offset);
    decTargetCode = targetCode(onset:offset);
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
    D(i).targetCode = decTargetCode;
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
initPertInd = find(decoderNum == perturbedDecoderNum, 1, 'first');
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

% Plot performance
if plotPerformance
    plotSessionPerformance(D);
end

% Determine whether or not to calculate performance statistics.  For all
% single-day experiments this should be done, but multi-day experiments
% don't have the same block design.
switch A.protocol
    case{'C_1.0','D_1.0','E_1.0'}
        calcStatistics = false;
        plotTrajectories = false;
end % Otherwise, use defaults (specified above)

if plotTrajectories
    % Plot trajectories
    trajBlockInds.baselineEval = baselineEvalInd;
    trajBlockInds.initPerturbation = initPertInd;
    trajBlockInds.stitchingEval = stitchingEvalInd;
    trajBlockInds.perturbationEval = perturbationEvalInd;
    plotNeilsenTrajectories(TDcell,trajBlockInds,E,saveFigs,saveLoc);
end

% Plot average success rate and acquisition time
if calcStatistics
    % Get indices for evaluation blocks and number of trials per eval block.
    % Also specify how many trials to discard at the beginning of each block --
    % this prevents learning effects from influencing the results
    evalInds = [baselineEvalInd stitchingEvalInd perturbationEvalInd baselineWashoutInd];
    blockOffset = [0 0 0 0]; % How many trials to discard at the beginning of the block.
    nEvalBlocks = length(evalInds);
    evalTrials = nan(nEvalBlocks,1);
    for i = 1:nEvalBlocks
        evalTrials(i) = D(evalInds(i)).nTrials - blockOffset(i);
    end
    
    % Baseline washout data typically lasted longer than 128 trials.  To
    % keep analysis consistent with other blocks and to limit the effect of
    % motivation, limit this block to 128 trials.  Also do the same for the
    % baseline evaluation block in case it went longer than 128 trials.
    maxEvalTrials = 128;
    evalTrials(end) = min(evalTrials(end),maxEvalTrials);
    
    blockTrials = floor(evalTrials/blockSize) * blockSize;      % Number of trials to analyze (only analyze blocks)
    successTrials = nan(nEvalBlocks,1);
    successRateMean = nan(nEvalBlocks,1);
    acquireTimeMean = nan(nEvalBlocks,1);
    acquireTimeStd = nan(nEvalBlocks,1);
    rewardRate = nan(nEvalBlocks,1);
    successRateData = cell(nEvalBlocks,1);
    acquireTimeDataSuccess = cell(nEvalBlocks,1);
    acquireTimeDataAll = cell(nEvalBlocks,1);
    targSR = cell(nEvalBlocks,1);
    targAT = cell(nEvalBlocks,1);
    targRR = cell(nEvalBlocks,1);
    targATQuant = cell(nEvalBlocks,1);
    blockNames = cell(nEvalBlocks,1);
    nTrialsTotal = nan(nEvalBlocks,1);
    nTrialsSuccess = nan(nEvalBlocks,1);

    % Get data for stitching block
    blockMask = ismember({D.blockNotes},'Stitching');
    if ~isempty(blockMask)
        % Get data and initialize matrices
        Dtemp = D(blockMask);
        nBlockTrials = length(Dtemp);
        updateNum = nan(1,nBlockTrials);
        nStitchingTrials = nan(1,nBlockTrials);     % Number of trials for each update
        sr_block = nan(1,nBlockTrials);             % Success rate for each block
        atAllMean_block = nan(1,nBlockTrials);      % Avg acquisition time for each block (all trials)
        atSuccessMean_block = nan(1,nBlockTrials);  % Avg acquisition time for each block (suc trials)
        tar_block = nan(1,nBlockTrials);            % Target acquisition rate
        
        % Loop over unique decoders and get data
        for i = 1:nBlockTrials
            % Get data for current block
            successCodeTemp = Dtemp(i).successCode';
            acquireTimeTemp = Dtemp(i).acquireTime';
            nStitchingTrials(i) = length(successCodeTemp);
            scMask = logical(successCodeTemp);
            
            % Put transposed data back into decoder structure array -- this
            % allows the success code and trial data to be easily
            % concatenated across all trials
            Dtemp(i).updateNum = ones(1,length(Dtemp(i).successCode)) * ...
                Dtemp(i).decoderNum;
            Dtemp(i).successCode = successCodeTemp;
            Dtemp(i).acquireTime = acquireTimeTemp;
            
            % Calculate stitching block performance metrics.  Compared to
            % normal 'block' metrics, stitching block metrics are
            % calculated for *each* stabilized decoder
            sr_block(i) = sum(successCodeTemp)/nStitchingTrials(i);
            atAllMean_block(i) = mean(acquireTimeTemp);
            atSuccessMean_block(i) = mean(acquireTimeTemp(scMask));
            tar_block(i) = sum(successCodeTemp)/sum(acquireTimeTemp);
        end
        % Get data for each trial in the stitching block
        stitchingUpdateNum = [Dtemp.updateNum];
        stitchingUpdateNum = stitchingUpdateNum - stitchingUpdateNum(1);
        successRateStitching = [Dtemp.successCode];
        acquireTimeStitching = [Dtemp.acquireTime];
        blockSuccessRateStitching = sr_block;
        switch acqTimeOption
            case 'successOnly'
                blockAcquireTimeStitching = atSuccessMean_block;
            case 'allTrials'
                blockAcquireTimeStitching = atAllMean_block;
        end
        blockTrialsStitching = nStitchingTrials;
        blockRewardRateStitching = tar_block;
    else
        successRateStitching = [];
        acquireTimeStitching = [];
        blockSuccessRateStitching = [];
        blockAcquireTimeStitching = [];
        blockRewardRateStitching = [];
    end
    
    % Loop over blocks to get data for evaluation blocks
    for i = 1:nEvalBlocks
        % Get success rate data, acquisition time data, and indices to analyze
        blockInd = evalInds(i);
        onset = blockOffset(i) + 1;
        
        % Updated 2017.11.17
        % Change the block data to use all trials for the evaluation blocks
        % instead of limiting things to blocks of 16 trials.  We should use
        % all available data when calculating performance metrics.
        
        %offset = onset + blockTrials(i) - 1;
        offset = onset + evalTrials(i) - 1; % Instead of 'blockTrials'
        
        successCodeTemp = D(blockInd).successCode(onset:offset);
        acquireTimeAllTemp = D(blockInd).acquireTime(onset:offset);

        % Calculate success rate
        successTrials(i) = sum(successCodeTemp);
        nTrialsTotal(i) = length(successCodeTemp);
        successRateData{i} = successCodeTemp;
        successRateMean(i) = mean(successCodeTemp);

        % Calculate reward rate
        rewardRate(i) = successTrials(i)/sum(acquireTimeAllTemp);
        
        % Calculate acquire time -- use successful trials ONLY
         % Select acquisition time data to calculate block average
        switch acqTimeOption
            case 'successOnly'
                acquireTimeTemp = acquireTimeAllTemp(logical(successCodeTemp));
            case 'allTrials'
                acquireTimeTemp = acquireTimeAllTemp;
            otherwise
                error('Incorrect acquisition time option specified.')
        end
        nTrialsSuccess(i) = sum(successCodeTemp);
        acquireTimeDataAll{i} = acquireTimeAllTemp;
        acquireTimeDataSuccess{i} = acquireTimeAllTemp(logical(successCodeTemp));
        acquireTimeMean(i) = mean(acquireTimeTemp);
        acquireTimeStd(i) = std(acquireTimeTemp);

        % Calculate success rate and acquisition time by target
        tC = D(blockInd).targetCode(onset:offset);
        uniTC = unique(tC);
        nTarg = length(uniTC);
        targSuccRate = nan(nTarg,1);        % Success rate
        targAcqTime = nan(nTarg,1);         % Mean acquisition time
        targRewardRate = nan(nTarg,1);
        q = [.25 .5 .75];                   % Quantiles to analyze
        targAcqTimeQuant = nan(nTarg,3);    % Acquisition time quantiles
        for j = 1:nTarg
            % Get mask for current target
            targMask = (tC == uniTC(j));
            aTMask = targMask(logical(successCodeTemp));
            tempSR = successCodeTemp(targMask);
            tempAT = acquireTimeAllTemp(targMask);
            targRewardRate(j) = sum(tempSR)/sum(tempAT);
            tempAT = acquireTimeTemp(aTMask);
            targSuccRate(j) = nanmean(tempSR);
            targAcqTime(j) = nanmean(tempAT);
            targAcqTimeQuant(j,:) = quantile(tempAT,q);
        end
        targSR{i} = targSuccRate;
        targAT{i} = targAcqTime;
        targRR{i} = targRewardRate;
        targATQuant{i} = targAcqTimeQuant;
        
        % Get block names
        blockNames{i} = D(blockInd).blockNotes;
    end

    % Calculate success rate statistics
    pValsSuccessRate = ones(nEvalBlocks);
    for i = 1:nEvalBlocks
        for j = (i+1):nEvalBlocks
            N1 = blockTrials(i); % Number of trials
            N2 = blockTrials(j);
            M1 = successTrials(i); % Number of successful trials
            M2 = successTrials(j);
            p = binomialTest(M1,N1,M2,N2);
            pValsSuccessRate(i,j) = p;
        end
    end

    % Calculate acquisition time statistics
    pValsAcquireTime = ones(nEvalBlocks);
    for i = 1:nEvalBlocks
        for j = (i+1):nEvalBlocks
            a1 = acquireTimeDataSuccess{i};
            a2 = acquireTimeDataSuccess{j};
            [~,p] = ttest2(a1,a2);
            pValsAcquireTime(i,j) = p;
        end
    end

    % Pack up data for evaluation blocks
    E.nEvalBlocks = nEvalBlocks;
    E.nTrialsTotal = nTrialsTotal;
    E.nTrialsSuccess = nTrialsSuccess;
    E.successRateData = successRateData;
    E.successRateMean = successRateMean;
    E.acquireTimeDataAll = acquireTimeDataAll;
    E.acquireTimeDataSuccess = acquireTimeDataSuccess;
    E.acquireTimeMean = acquireTimeMean;
    E.acquireTimeStd = acquireTimeStd;
    E.rewardRate = rewardRate;
    E.targetSuccessRate = targSR;
    E.targetAcquireTime = targAT;
    E.targetRewardRate = targRR;
    E.targetAcquireTimeQuant = targATQuant;
    E.pValsSuccessRate = pValsSuccessRate;
    E.pValsAcquireTime = pValsAcquireTime;
    E.blockNames = blockNames;
    
    % Pack up data for stitching
    E.stitchingUpdateNum = stitchingUpdateNum;
    E.successRateStitching = successRateStitching;
    E.acquireTimeStitching = acquireTimeStitching;
    E.blockTrialsStitching = blockTrialsStitching;
    E.blockSuccessRateStitching = blockSuccessRateStitching;
    E.blockAcquireTimeStitching = blockAcquireTimeStitching;
    E.blockRewardRateStitching = blockRewardRateStitching;

    if plotStatistics
        plotStitchingStatistics(E,saveFigs,saveLoc)
    end

    % Determine Perturbation type
    isOffset = E.pertParams.nOffsetChs > 0;
    isSilent = E.pertParams.nSilentChs > 0;
    isSwap = E.pertParams.nSwapChs > 0;

    if isOffset && ~(isSilent || isSwap)
        pertType = 'offset';
    elseif isSilent && ~(isOffset || isSwap)
        pertType = 'silence';
    elseif isSwap && ~(isSilent || isOffset)
        pertType = 'swap';
    elseif ~(isOffset || isSilent || isSwap)
        pertType = 'control';
    else
        pertType = 'multiple';
    end
    E.pertType = pertType;

    % Save experiment meta data
    if saveData
        saveLoc = [E.dataLocBase '/ExperimentData/' E.subject E.dataset ...
            '_ExperimentData.mat'];
        save(saveLoc,'E')
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
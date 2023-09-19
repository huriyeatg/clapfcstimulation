% generateStimuliWithLEDMaskers
% This function generates .mat files for all atimulation vectors for Packer
% I/O and saves them to the folder called finalVersion.
% Three files will be generated;
%               1) Stimulation (30 trials for stimulation, 6 no stimulation)
%               2) LED maskers (36 trials, all trials + 20% control
%               trials)
%               3) Shutter (36 trials, starts 50ms before the stimulation,
%               as takes time to close)
%               trial,
% There will be few minutes spontaneous recording to start with.

tic
%% parameters
params.savePath = fullfile(root_path,'onePhotonStimulation','stimuli_onePhotonStimulation',['StimuliForTraining_',date]);
if ~exist(params.savePath,'dir')
    mkdir(params.savePath)
end
params.fs = 20000;% as in the Packer I/0
params.nTrialsCueReward    = 50;
params.nTrialsCueNoReward  = 10;
params.cueDuration = 2; %sec
params.rewardDuration = 0.05; %sec
params.amp           = 2.45;  % Hz
params.cue_amp      = params.amp*2;
params.reward_amp      = params.amp*2;
params.totalTrialDuration = params.cueDuration + params.rewardDuration;
params.itiCutOff = 12; %sec
params.itiStart  = 20;
fs = params.fs;

%% LED cue signal
t = 1/fs:(1/fs):params.cueDuration;
trial_cueSignal = (ones(size(t))*params.cue_amp);%+rand(size(t))/10;
% last timepoint come back to zero
trial_cueSignal = [trial_cueSignal(1:end-1),zeros(1,(params.rewardDuration*fs +1 ))];

%% Water/reward delivery signal
t = 1/fs:(1/fs):params.rewardDuration;
trial_rewardSignal = (ones(size(t))*params.reward_amp);%+rand(size(t))/10;
% last timepoint come back to zero
trial_rewardSignal = [zeros(1,(params.cueDuration*fs)),trial_rewardSignal(1:end-1),0];

%% plot the trials
t = 1/fs:(1/fs):params.totalTrialDuration;
figure; 
subplot(3,1,1)
plot(t,trial_cueSignal); title('Cue Signal');box off; ylim([0 5])
ylabel('Volt')
subplot(3,1,2)
plot(t,trial_rewardSignal); title('Reward Signal');box off; ylim([0 5])
ylabel('Volt')
xlabel('Time (sec)')

saveas(gcf,fullfile( params.savePath, ['trialFig',date,'.fig']), 'fig');

params.trial_cueSignal    = trial_cueSignal;
params.trial_rewardSignal    = trial_rewardSignal;
%% %%%%%% Create the block
% Add spontaneous recording

params.itiTimes = params.itiStart + (params.itiCutOff-params.itiStart) .* rand(params.nTrialsCueReward + params.nTrialsCueNoReward,1);

params.trialIDfornoRewardtrials = randi(params.nTrialsCueReward + params.nTrialsCueNoReward,  params.nTrialsCueNoReward,1);
% Create 120 trials

blockCueSignal = [];
blockRewardSignal = [];
for k=1:(params.nTrialsCueReward + params.nTrialsCueNoReward)
    % Create no-signal vector for iti
    temp = params.itiTimes(k,1);
    t = 1/fs:(1/fs):temp;
    itiSignal = zeros(size(t));%+rand(size(t))/10;   
    
    if ismember(k, params.trialIDfornoRewardtrials) % means no reward 
        t = 1/fs:(1/fs):(params.rewardDuration + params.cueDuration);
        temp_rewardSignal = zeros(size(t));%+rand(size(t))/10;
    else
        temp_rewardSignal = trial_rewardSignal;
    end
        
    blockCueSignal     = [ blockCueSignal, trial_cueSignal, itiSignal];
    blockRewardSignal  = [ blockRewardSignal, temp_rewardSignal, itiSignal];   
end

params.blockCueSignal     = blockCueSignal;
params.blockRewardSignal  = blockRewardSignal;

params.totalBlockTime = size(blockCueSignal,2)/fs; % sec
%% %%%%%% Plot the block
t = 1/fs:(1/fs):params.totalBlockTime;
figure; 
subplot(3,1,1)
plot(t,blockCueSignal); title('Cue Signal');box off; ylim([0 5])
ylabel('Volt')
subplot(3,1,2)
plot(t,blockRewardSignal); title('Reward Signal');box off; ylim([0 5])
ylabel('Volt')
xlabel('Time (sec)')
saveas(gcf,fullfile( params.savePath, ['blockFig',date,'.fig']), 'fig');
print(gcf,fullfile( params.savePath, ['blockFig',date,'.png']), '-dpng');

% save the timeSeries in .dat format
disp('Files are saving, might take few seconds')
fid=fopen(fullfile( params.savePath, ['cueSignal_',date,'.dat']),'w','l');
fwrite(fid,blockCueSignal,'double');fclose(fid);

fid=fopen(fullfile( params.savePath, ['rewardSignal_',date,'.dat']),'w','l');
fwrite(fid,blockRewardSignal,'double');fclose(fid);

save( fullfile( params.savePath, ['infoStim',date,'.mat']), 'params')
toc
fprintf('Files are saved in %s\n',params.savePath)
close all
% 
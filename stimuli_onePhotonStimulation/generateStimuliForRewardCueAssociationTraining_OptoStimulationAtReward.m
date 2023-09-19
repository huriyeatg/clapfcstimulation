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
params.savePath = fullfile(root_path,'onePhotonStimulation','stimuli_onePhotonStimulation',['StimuliForTraining_OptoAtReward',date]);
if ~exist(params.savePath,'dir')
    mkdir(params.savePath)
end
params.fs = 20000;% as in the Packer I/0
params.nTrialsCueReward      = 60;
params.nTrialsCueRewardStim  = 30;
params.nTrialsCueNoReward    = 10;
params.cueDuration = 2; %sec
params.pulseDuration = params.cueDuration;
params.rewardDuration = 0.05; %sec
params.amp           = 2.45;  % Hz
params.cue_amp      = params.amp*2;
params.reward_amp      = params.amp*2;
params.totalTrialDuration = params.cueDuration + params.pulseDuration;
params.itiCutOff = 12; %sec
params.itiStart  = 20;
fs = params.fs;
params.offset        = 2.45;
params.amp           = 2.45;  % Hz
params.dutyCycleVal  = 50; % Percentage
params.freq          = 40; % Hz


%% stimulation signal
%Generate a 25 Hz square wave sampled at 1 kHz for 70 ms.
%Specify a duty cycle of 5%. Add white Gaussian noise with a variance of
%1/100.
t = 1/fs:(1/fs):params.pulseDuration;
trial_stimSignal = params.offset + (params.amp*square(2*pi*params.freq*t,params.dutyCycleVal));%+rand(size(t))/10;
% last timepoint has another pulse to start - clean it!
trial_stimSignal = [zeros(1,(params.cueDuration*fs)),trial_stimSignal(1:end-1),0];

%% LED cue signal
t = 1/fs:(1/fs):params.cueDuration;
trial_cueSignal = (ones(size(t))*params.cue_amp);%+rand(size(t))/10;
% last timepoint come back to zero
trial_cueSignal = [trial_cueSignal(1:end-1),zeros(1,(params.pulseDuration*fs +1 ))];

%% Water/reward delivery signal
t = 1/fs:(1/fs):params.rewardDuration;
trial_rewardSignal = (ones(size(t))*params.reward_amp);%+rand(size(t))/10;
% last timepoint come back to zero
trial_rewardSignal = [zeros(1,(params.cueDuration*fs)),trial_rewardSignal(1:end-1),zeros(1,((params.pulseDuration-params.rewardDuration)*fs +1 ))];

%% plot the trials
t = 1/fs:(1/fs):params.totalTrialDuration;
figure; 
subplot(3,1,1)
plot(t,trial_cueSignal); title('Cue Signal');box off; ylim([0 5])
ylabel('Volt')
subplot(3,1,2)
plot(t,trial_stimSignal); title('Photostimulation Signal');box off; ylim([0 5])
ylabel('Volt')
subplot(3,1,3)
plot(t,trial_rewardSignal); title('Reward Signal');box off; ylim([0 5])
ylabel('Volt')
xlabel('Time (sec)')

saveas(gcf,fullfile( params.savePath, ['trialFig',date,'.fig']), 'fig');

params.trial_cueSignal    = trial_cueSignal;
params.trial_rewardSignal    = trial_rewardSignal;
params.trial_timSignal   = trial_stimSignal;
%% %%%%%% Create the block
% Add spontaneous recording

params.itiTimes = params.itiStart + (params.itiCutOff-params.itiStart) .* rand(params.nTrialsCueReward + params.nTrialsCueRewardStim,1);

params.trialIDforStimTrials = randi(params.nTrialsCueReward + params.nTrialsCueRewardStim,  params.nTrialsCueRewardStim,1);

% Create block with all trials
blockCueSignal = [];
blockRewardSignal = [];
blockStimulationSignal = [];
for k=1:(params.nTrialsCueReward + params.nTrialsCueRewardStim)
    % Create no-signal vector for iti
    temp = params.itiTimes(k,1);
    t = 1/fs:(1/fs):temp;
    itiSignal = zeros(size(t));%+rand(size(t))/10;   
    
    if ismember(k, params.trialIDforStimTrials) % means no stimulation 
        t = 1/fs:(1/fs):params.totalTrialDuration;
        temp_stimSignal = zeros(size(t));%+rand(size(t))/10;
    else
        temp_stimSignal = trial_stimSignal;
    end
        
    blockCueSignal     = [ blockCueSignal, trial_cueSignal, itiSignal];
    blockRewardSignal  = [ blockRewardSignal, trial_rewardSignal, itiSignal];   
    blockStimulationSignal  = [ blockStimulationSignal, temp_stimSignal, itiSignal]; 
end

params.blockCueSignal     = blockCueSignal;
params.blockRewardSignal  = blockRewardSignal;
params.blockStimulationSignal  = blockStimulationSignal;

params.totalBlockTime = size(blockCueSignal,2)/fs; % sec
%% %%%%%% Plot the block
t = 1/fs:(1/fs):params.totalBlockTime;
figure; 
subplot(3,1,1)
plot(t,blockCueSignal); title('Cue Signal');box off; ylim([0 5])
ylabel('Volt')
subplot(3,1,2)
plot(t,blockStimulationSignal); title('PhotoStimulation Signal');box off; ylim([0 5])
ylabel('Volt')
subplot(3,1,3)
plot(t,blockRewardSignal); title('Reward Signal');box off; ylim([0 5])
ylabel('Volt')
xlabel('Time (sec)')
saveas(gcf,fullfile( params.savePath, ['blockFig',date,'.fig']), 'fig');
print(gcf,fullfile( params.savePath, ['blockFig',date,'.png']), '-dpng');

%% save the timeSeries in .dat format
disp('Files are saving, might take few seconds')
fid=fopen(fullfile( params.savePath, ['cueSignal_',date,'.dat']),'w','l');
fwrite(fid,blockCueSignal,'double');fclose(fid);

fid=fopen(fullfile( params.savePath, ['rewardSignal_',date,'.dat']),'w','l');
fwrite(fid,blockRewardSignal,'double');fclose(fid);

fid=fopen(fullfile( params.savePath, ['stimulationSignal_',date,'.dat']),'w','l');
fwrite(fid,blockStimulationSignal,'double');fclose(fid);

save( fullfile( params.savePath, ['infoStim',date,'.mat']), 'params')
toc
fprintf('Files are saved in %s\n',params.savePath)
close all
% 
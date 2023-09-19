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
params.savePath = fullfile(root_path,'onePhotonStimulation','stimuli_onePhotonStimulation',['Stimuli_',date]);
if ~exist(params.savePath,'dir')
    mkdir(params.savePath)
end
params.fs = 20000;% as in the Packer I/0
params.nTrialsNoLight = 20;
params.nTrialsLight   = 30;
params.sponDur       = 180;  % sec
params.shutterDelay  = 0.05; % 50ms 
params.pulseDuration = 0.25;  % in sec, 0.02 = 20ms
params.offset        = 2.45;
params.amp           = 2.45;  % Hz
params.dutyCycleVal  = 50; % Percentage
params.freq          = 40; % Hz
params.masker_amp      = params.amp*2;
params.maskerDuration  = params.pulseDuration;
params.shutter_amp     = params.amp*2;
params.itiCutOff = 9; %sec
params.itiStart  = 15;

fs = params.fs;
%% spontaneous recording signal
t = 1/fs:(1/fs):params.sponDur;
sponSignal = zeros(size(t));%+rand(size(t))/10

%% stimulation signal
%Generate a 25 Hz square wave sampled at 1 kHz for 70 ms.
%Specify a duty cycle of 5%. Add white Gaussian noise with a variance of
%1/100.
t = 1/fs:(1/fs):params.pulseDuration;
trial_stimSignal = params.offset + (params.amp*square(2*pi*params.freq*t,params.dutyCycleVal));%+rand(size(t))/10;
% last timepoint has another pulse to start - clean it!
trial_stimSignal = [trial_stimSignal(1:end-1),0];
%dutycycle(trial_stimSignal,fs);


%% LED masker signal
t = 1/fs:(1/fs):params.maskerDuration;
trial_maskerSignal = (ones(size(t))*params.masker_amp);%+rand(size(t))/10;
% last timepoint come back to zero
trial_maskerSignal = [trial_maskerSignal(1:end-1),0];

%% shutter signal
params.shutterDuration = params.pulseDuration + params.shutterDelay; 
t = 1/fs:(1/fs):params.shutterDuration;
trial_shutterSignal = (ones(size(t))*params.shutter_amp);%+rand(size(t))/10;
% last timepoint come back to zero
trial_shutterSignal = [trial_shutterSignal(1:end-1),0];

%!! UPDATE stimulation and LED signal& add shutter delay before the 
%stimulation and LED starts
t = 1/fs:(1/fs):params.shutterDelay;
delayTime = zeros(size(t));%+rand(size(t))/10;
trial_stimSignal   = [delayTime,trial_stimSignal];
trial_maskerSignal = [delayTime,trial_maskerSignal];

%% plot the trials
t = 1/fs:(1/fs):params.shutterDuration;
figure; 
subplot(3,1,1)
plot(t,trial_stimSignal); title('Stimulation pulses');box off; ylim([0 5])
ylabel('Volt')
subplot(3,1,2)
plot(t,trial_shutterSignal); title('Shutter Pulse');box off; ylim([0 5])
ylabel('Volt')
subplot(3,1,3)
plot(t,trial_maskerSignal); title('Masker LED light');box off; ylim([0 5])
ylabel('Volt')
xlabel('Time (sec)')

saveas(gcf,fullfile( params.savePath, ['trialFig',date,'.fig']), 'fig');

params.trial_stimSignal    = trial_stimSignal;
params.trial_shutterSignal = trial_shutterSignal;
params.trial_maskerSignal  = trial_maskerSignal;
params.sponSignal          = sponSignal;

%% %%%%%% Create the block
% Add spontaneous recording
stimulationSignal = [sponSignal];
shutterSignal = [sponSignal];
maskerSignal = [sponSignal];

params.itiTimes = params.itiStart + (params.itiCutOff-params.itiStart) .* rand(params.nTrialsLight+params.nTrialsNoLight,1);

% Create 20 trials without light
eventTimes =[];eventID =[];
lightID = 0;   stimulationID = 1;
for k=1: params.nTrialsNoLight
    temp = params.itiTimes(k,1);
    t = 1/fs:(1/fs):temp;
    itiSignal = zeros(size(t));%+rand(size(t))/10;

    tempstimSignal =  trial_stimSignal;
    
    t = 1/fs:(1/fs):params.shutterDuration;
    tempmaskerSignal = zeros(size(t));%+rand(size(t))/10;
       

    eventTimes       = [ eventTimes, size(stimulationSignal,2)/fs];
    eventID          = [eventID; lightID, stimulationID]; % light, stimulation
    stimulationSignal= [ stimulationSignal, tempstimSignal, itiSignal];
    shutterSignal    = [ shutterSignal, trial_shutterSignal, itiSignal];
    maskerSignal     = [ maskerSignal, tempmaskerSignal, itiSignal];
    
end
sum(eventID)

%Create 36 trials with light only 20 percent of them without light
tempmaskerSignal = trial_maskerSignal; % all trials will have light
tInd =0;   lightID = 1;
for k=1: params.nTrialsLight
    temp = params.itiTimes(k+params.nTrialsNoLight,1);
    t = 1/fs:(1/fs):temp;
    itiSignal = zeros(size(t));%+rand(size(t))/10;
    % 20% for no stimulation
    if rand()<=0.33
        tInd = tInd+1;
        stimulationID = 0;% no stimution
        t = 1/fs:(1/fs):params.shutterDuration;
        tempstimSignal = zeros(size(t));%+rand(size(t))/10;
    else
        stimulationID = 1;%stimulation
        tempstimSignal =  trial_stimSignal;
    end

    eventTimes       = [ eventTimes, size(stimulationSignal,2)/fs];
    eventID          = [eventID; lightID, stimulationID]; % light, stimulation
    stimulationSignal= [ stimulationSignal, tempstimSignal, itiSignal];
    shutterSignal    = [ shutterSignal, trial_shutterSignal, itiSignal];
    maskerSignal     = [ maskerSignal, tempmaskerSignal, itiSignal];
    
end
tInd
sum(eventID)
params.totalBlockTime = size(stimulationSignal,2)/fs;
params.eventTimes     = eventTimes;
params.eventID        = eventID;

% add eventTriggerSignal - first point of all events ( stimulation + light)
eventTriggerSignal    = zeros(1,size(maskerSignal,2));
[pks, loc] = findpeaks(maskerSignal+stimulationSignal,'MinPeakDistance',size(trial_stimSignal,2));
eventTriggerSignal(loc) = params.masker_amp;

%% %%%%%% Plot the block
t = 1/fs:(1/fs):params.totalBlockTime;
figure; 
subplot(3,1,1)
plot(t,stimulationSignal); title('Stimulation pulses');box off; ylim([0 5])
ylabel('Volt')
subplot(3,1,2)
plot(t,shutterSignal); title('Shutter Pulse');box off; ylim([0 5])
ylabel('Volt')
subplot(3,1,3)
plot(t,maskerSignal); title('Masker LED light');box off; ylim([0 5])
ylabel('Volt')
xlabel('Time (sec)')
saveas(gcf,fullfile( params.savePath, ['blockFig',date,'.fig']), 'fig');
print(gcf,fullfile( params.savePath, ['blockFig',date,'.png']), '-dpng');

% save the timeSeries in .dat format
disp('Files are saving, might take few seconds')
fid=fopen(fullfile( params.savePath, ['stimulationSignal_',date,'.dat']),'w','l');
fwrite(fid,stimulationSignal,'double');fclose(fid);

fid=fopen(fullfile( params.savePath, ['shutterSignal_',date,'.dat']),'w','l');
fwrite(fid,shutterSignal,'double');fclose(fid);

fid=fopen(fullfile( params.savePath, ['maskerSignal_',date,'.dat']),'w','l');
fwrite(fid,maskerSignal,'double');fclose(fid);

fid=fopen(fullfile( params.savePath, ['eventTriggerSignal',date,'.dat']),'w','l');
fwrite(fid,eventTriggerSignal,'double');fclose(fid);

save( fullfile( params.savePath, ['infoStim',date,'.mat']), 'params','stimulationSignal','shutterSignal','maskerSignal','eventTriggerSignal')
toc
fprintf('Files are saved in %s\n',params.savePath)
close all
% 
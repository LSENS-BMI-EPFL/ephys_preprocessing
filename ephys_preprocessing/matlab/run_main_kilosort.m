%% you need to change most of the paths in this block
function [] = run_main_kilosort(bin_data_path, fs, temp_data_path)

addpath(genpath('C:\Users\bisi\Kilosort\Kilosort-2.0')) % path to kilosort folder
addpath('C:\Users\bisi\Github\npy-matlab') % for converting to Phy

rootZ = bin_data_path; % raw data binary file is in this folder
%rootH = rootZ;  % path to temporary binary file (same size as data, should be on fast SSD)
outputFolder = [rootZ, '\kilosort2'];

pathToYourConfigFile = 'C:\Users\bisi\Kilosort\Kilosort-2.0\configFiles'; % take from Github folder and put it somewhere else (together with the master_file)
chanMapFile = 'neuropixPhase3B1_kilosortChanMap.mat'; % check

run(fullfile(pathToYourConfigFile, 'configFile384.m')) % this sets all parameters

ops.trange = [0 Inf]; % time range to sort
ops.NchanTOT    = 385; % total number of channels in your recording
ops.fs = double(fs); % set sampling rate manually
ops.fproc       = fullfile(outputFolder, 'temp_wh.dat'); % proc file on a fast SSD
ops.chanMap = fullfile(pathToYourConfigFile, chanMapFile);


%% this block runs all the steps of the algorithm
fprintf('Looking for data inside %s \n', rootZ)

% is there a channel map file in this folder?
fs = dir(fullfile(rootZ, 'chan*.mat'));
if ~isempty(fs)
    ops.chanMap = fullfile(rootZ, fs(1).name);
end

% find the binary file
%fs = [dir(fullfile(rootZ, 'corrected*.bin')) dir(fullfile(rootZ, '*.dat'))]; .dat file are preprocessed data
fs = [dir(fullfile(rootZ, '*corrected*.bin'))];
ops.fbinary = fullfile(rootZ, fs(1).name);
disp(['Binary file to spike-sort: ' ops.fbinary]);

% preprocess data to create temp_wh.dat
rez = preprocessDataSub(ops);

% time-reordering as a function of drift
rez = clusterSingleBatches(rez);

% saving here is a good idea, because the rest can be resumed after loading rez
save(fullfile(outputFolder, 'rez.mat'), 'rez', '-v7.3');

% main tracking and template matching algorithm
rez = learnAndSolve8b(rez);

% OPTIONAL: remove double-counted spikes - solves issue in which individual spikes are assigned to multiple templates.
% See issue 29: https://github.com/MouseLand/Kilosort2/issues/29
rez = remove_ks2_duplicate_spikes(rez);

% final merges
rez = find_merges(rez, 1);

% final splits by SVD
rez = splitAllClusters(rez, 1);

% final splits by amplitudes
rez = splitAllClusters(rez, 0);

% decide on cutoff
rez = set_cutoff(rez);

fprintf('found %d good units \n', sum(rez.good>0))

% write to Phy

fprintf('Saving results to Phy  \n')
rezToPhy(rez, outputFolder);


end



%% Pipeline
% creates a table listing stimuli pairs and corresponding reponses (2 =
% upper/first pair is more similar, 3 = lower/second pair is more similar)
% for ALL days per subject

clear all
close all

clc


%% SET PATH
%
mainpath = 'C:\Users\rumpf\Documents\perceptual_space';
%

datapath = [mainpath filesep 'logfiles' filesep 'data'];

% Set working dictonary
cd(mainpath);
addpath ('C:\Users\rumpf\Documents\perceptual_space\maximum likelihood difference scaling');

%% DEFINE PARAMETERS

%Subject_Index = {'02', '03'};, '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15'};
Subject_Index = {'15'};
%Subject_Index = {'09'; '10'; '11'; '12'; '13'; '14'; '15'};

for i = length(Subject_Index)
    load([datapath filesep 'ps_' Subject_Index{i} '_day1_stim_resp.mat'], 'stim_resp');
    day1 = stim_resp;

    load([datapath filesep 'ps_' Subject_Index{i} '_day2_stim_resp.mat'], 'stim_resp');
    day2 = stim_resp;

    stim_resp_all = [day1;day2];

    save([datapath filesep 'ps_' Subject_Index{i} '_alldays_stim_resp.mat'], 'stim_resp_all');
end
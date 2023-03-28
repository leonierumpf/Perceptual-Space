function mean_space = align_subjects
%
%Schritte
% 1. Alle perceptual spaces laden (jeweils eine 8x2 matrix)
% 2. Optimalen Kreis als Positionen definieren
% 3. Alle perceptual spaces an optimalen Kreis procrustes alignen
% 4. Mean über die alignten perceptual spaces berechnen
clear all
close all

clc


%% SET PATH
%
mainpath = 'C:\Users\rumpf\Documents\perceptual_space';
%
%stimpath = [mainpath filesep 'logfiles' filesep 'stimlists'];

%datapath = [mainpath filesep 'logfiles' filesep 'data'];

% Set working dictonary
cd(mainpath);
addpath ('C:\Users\rumpf\Documents\perceptual_space\logfiles\data\best_positions_day2')

%% 1. alle Spaces von allen subjects laden
%matrix mit dimensionen: Stimuli x Dimensionen x Subjects
subject_spaces = randn(8, 2, 16); % insert best position of each participant 

path_directory = 'C:\Users\rumpf\Documents\perceptual_space\logfiles\data\best_positions_day2';
original_files = dir([path_directory filesep '*.mat']);
bp = randn(8,2);

for k =1:length(original_files)
    
    filename = [path_directory '/' original_files(k).name];
    load(original_files(k).name, '-mat')
    bp(:,1) = best_position(1,:).';
    bp(:,2) = best_position(2,:).';
    subject_spaces(:,:,k) = bp;

end 


%% 2. Optimalen Kreis definieren
% Idee: Wie sieht der Space aus, wenn Probanden nur die Gradzahl verwenden,
% um die Ähnlichkeit von Stimuli zu definieren
% Winkel der Stimuli in Radians
angles = deg2rad(linspace(0, 315, 8));
space = [cos(angles); sin(angles)]';

%% 3. Alle subject spaces an den Kreis alignen
aligned_spaces = zeros(size(subject_spaces));
for sub = 1:size(aligned_spaces, 3)
    [~, aligned_spaces(:, :, sub)] = procrustes(space, subject_spaces(:, :, sub));
end

%% 4. Average über die aligned spaces
mean_space = mean(aligned_spaces, 3);

% wenn mean_space wie ein Kreis aussieht, kannst du für die Auswertung
% einfach davon ausgehen, dass der perzeptuelle Raum in etwa ein Kreis ist.

scatter(mean_space(:, 1), mean_space(:, 2))
axis equal
for i = 1:8
        text(mean_space(i, 1), mean_space(i, 2), num2str(i));
    end
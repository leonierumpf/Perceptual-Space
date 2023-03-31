% This is the MLE function fitted to data from PS (stimuli: Gabor Patches).
% You need to insert your stim(ST)-resp(R) data, the number of dimensions you want
% to model your data in (ModelDimen) and how often this process should be repeated
% (n_iterations).
%
% Be mindful: This is an analysis of the PS per subject, for a pooled
% analyis refer to "align_subjects"
%
% ESTIMATE contains the coordinates of the perceptual scale values in MODELDIMENx8 
% matrix.
%
% Depends on MLDS_TriadLikelihood and MLDS_MLE

clear all
close all

clc


%% SET PATH
%
mainpath = 'C:\Users\rumpf\Documents\perceptual_space';
%
% in stimpath all the nice tables of each subjects are stored in a 216x5
% matrix
stimpath = 'C:\Users\rumpf\Documents\perceptual_space\logfiles\data\stim_resp';

savepath = 'C:\Users\rumpf\Documents\perceptual_space\best_positions';

% Set working dictonary
cd(mainpath);
% all function scripts are stored here:
addpath ('C:\Users\rumpf\Documents\perceptual_space\MLDS scripts');

%% DEFINE PARAMETERS D1
% these are subject to change
Subject_Index = { '02', '03', '04', '05', '06', '07', '08', '10', '11', '14', '15'};
ModelDimen = 1;
n_iterations = 100;

%% start loop day1

for si = 1:length(Subject_Index)

    stim_resp = load([stimpath filesep 'ps_' Subject_Index{si} '_day1_stim_resp.mat'], 'stim_resp'); % read in stimlists

    stim_resp = struct2array(stim_resp);
    
    ST = stim_resp(:,1:4);
    R = stim_resp(:,5);

    % run 100 iterations with random start values
    % if current one is better (i.e. lower negative log likelihood), keep
    % it
    % otherwise keep the old best one

    best_neg_ll = inf;
    best_position = nan;
    n_problems = 0;

    if ModelDimen == 1
        %Falls du später rumexperimentieren möchtest, parfor-loop über Subjects
        %um zu parallelisieren. Dafür brauchst du die parpool function
        for iter = 1:n_iterations
            fprintf('.')
            [Estimate,ExitFlag,Likelihood]=MLDS_MLE(ST,R,ModelDimen);
            if ExitFlag ~= 1;
                n_problems = n_problems + 1;
                continue
            end
            if Likelihood < best_neg_ll;
                best_neg_ll = Likelihood;
                best_position = reshape(Estimate(1:end-1), 1, 8); %Estimate(1:end-1) because we fixate two points

            end
        end

        fprintf('Did not converge in %d iterations  (%.1f %%)\n', n_problems, n_problems / n_iterations * 100);
        fprintf('Best average likelihood: %.2f %%\n', exp(-best_neg_ll) * 100);
        % fprintf('Best average likelihood: %.2f %%\n', best_neg_ll * 100);

        figure;
       % scatter(best_position(1, :), best_position(1, :));
         scatter(best_position(1, :), 0);
        for i = 1:8
           % text(best_position(1, i), best_position(1, i), num2str(i));
            text(best_position(1, i), 0, num2str(i));
        end

        savefig([savepath filesep '1D_best_positions_day1' filesep 'ps_' Subject_Index{si} '_1D_best_position_day1.fig']);
        save([savepath filesep '1D_best_positions_day1' filesep 'ps_' Subject_Index{si} '_1D_best_position_day1.mat']);
    
    elseif ModelDimen == 2

        %Falls du später rumexperimentieren möchtest, parfor-loop über Subjects
        %um zu parallelisieren. Dafür brauchst du die parpool function
        for iter = 1:n_iterations
            fprintf('.')
            [Estimate,ExitFlag,Likelihood]=MLDS_MLE(ST,R,ModelDimen);
            if ExitFlag ~= 1;
                n_problems = n_problems + 1;
                continue
            end
            if Likelihood < best_neg_ll;
                best_neg_ll = Likelihood;
                best_position = reshape(Estimate(1:end-1), 2, 8); %Estimate(1:end-1) because we fixate two points

            end
        end

        fprintf('Did not converge in %d iterations  (%.1f %%)\n', n_problems, n_problems / n_iterations * 100);
        fprintf('Best average likelihood: %.2f %%\n', exp(-best_neg_ll) * 100);
        % fprintf('Best average likelihood: %.2f %%\n', best_neg_ll * 100);

        figure;
        scatter(best_position(1, :), best_position(2, :));
        for i = 1:8
            text(best_position(1, i), best_position(2, i), num2str(i));
        end


        save([savepath filesep '2D_best_positions_day1' filesep 'ps_' Subject_Index{si} '2D_best_position_day1.mat'], '2D_best_position_d1');

    end
end


%% DEFINE PARAMETERS D2
% these are subject to change
Subject_Index = { '02', '03', '04', '05', '07', '09', '10', '11', '12', '14'};
ModelDimen = 2;
n_iterations = 50;

%% start loop day2
for si = 1:length(Subject_Index)

    % read in stimlists
    stim_resp = load([stimpath filesep 'ps_' Subject_Index{si} '_day2_stim_resp.mat'], 'stim_resp');
    
    stim_resp = struct2array(stim_resp);

    ST = stim_resp(:,1:4);
    R = stim_resp(:,5);

    % run 100 iterations with random start values
    % if current one is better (i.e. lower negative log likelihood), keep
    % it
    % otherwise keep the old best one

    best_neg_ll = inf;
    best_position = nan;
    n_problems = 0;

    if ModelDimen == 1
        %Falls du später rumexperimentieren möchtest, parfor-loop über Subjects
        %um zu parallelisieren. Dafür brauchst du die parpool function
        for iter = 1:n_iterations
            fprintf('.')
            [Estimate,ExitFlag,Likelihood]=MLDS_MLE(ST,R,ModelDimen);
            if ExitFlag ~= 1;
                n_problems = n_problems + 1;
                continue
            end
            if Likelihood < best_neg_ll;
                best_neg_ll = Likelihood;
                best_position = reshape(Estimate(1:end-1), 1, 8); %Estimate(1:end-1) because we fixate two points

            end
        end

        fprintf('Did not converge in %d iterations  (%.1f %%)\n', n_problems, n_problems / n_iterations * 100);
        fprintf('Best average likelihood: %.2f %%\n', exp(-best_neg_ll) * 100);
        % fprintf('Best average likelihood: %.2f %%\n', best_neg_ll * 100);

         figure;
       % scatter(best_position(1, :), best_position(1, :));
         scatter(best_position(1, :), 0);
        for i = 1:8
           % text(best_position(1, i), best_position(1, i), num2str(i));
            text(best_position(1, i), 0, num2str(i));
        end

        savefig([savepath filesep '1D_best_positions_day2' filesep 'ps_' Subject_Index{si} '1D_best_position_day2.fig']);
        save([savepath filesep '1D_best_positions_day2' filesep 'ps_' Subject_Index{si} '1D_best_position_day2.mat'], '1D_best_position_d2');

    elseif ModelDimen == 2

        %Falls du später rumexperimentieren möchtest, parfor-loop über Subjects
        %um zu parallelisieren. Dafür brauchst du die parpool function
        for iter = 1:n_iterations
            fprintf('.')
            [Estimate,ExitFlag,Likelihood]=MLDS_MLE(ST,R,ModelDimen);
            if ExitFlag ~= 1;
                n_problems = n_problems + 1;
                continue
            end
            if Likelihood < best_neg_ll;
                best_neg_ll = Likelihood;
                best_position = reshape(Estimate(1:end-1), 2, 8); %Estimate(1:end-1) because we fixate two points

            end
        end

        fprintf('Did not converge in %d iterations  (%.1f %%)\n', n_problems, n_problems / n_iterations * 100);
        fprintf('Best average likelihood: %.2f %%\n', exp(-best_neg_ll) * 100);
        % fprintf('Best average likelihood: %.2f %%\n', best_neg_ll * 100);

        figure;
        scatter(best_position(1, :), best_position(2, :));
        for i = 1:8
            text(best_position(1, i), best_position(2, i), num2str(i));
        end


        save([savepath filesep '2D_best_positions_day2' filesep 'ps_' Subject_Index{si} '2D_best_position_day2.mat'], '2D_best_position_d2');

    end
end
    save([datapath filesep 'ps_' Subject_Index{si} '_best_position_day1.mat'], 'best_position');
    save([datapath filesep 'ps_' Subject_Index{si} '_best_position_day2.mat'], 'best_position');
       %To do: Spaces zwischen Subjects alignen
       %procrustes (MATLAB funktion)
       %Für average über Subjects. Procrustes macht spaces so ähnlich wie
       %möglich, ohne die relative position von Stimuli zu verändern
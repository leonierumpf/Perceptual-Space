%% Pipeline
% creates a table listing stimuli pairs and corresponding reponses (2 =
% upper/first pair is more similar, 3 = lower/second pair is more similar)
% before you import data logfile (incl. responses) put them into .xlsx
% format

clear all
close all

clc


%% SET PATH
%
mainpath = 'C:\Users\rumpf\Documents\perceptual_space';
%
stimpath = [mainpath filesep 'logfiles' filesep 'stimlists']

datapath = [mainpath filesep 'logfiles' filesep 'data'];

% Set working dictonary
cd(mainpath);
addpath ('C:\Users\rumpf\Documents\perceptual_space\maximum likelihood difference scaling')

%% DEFINE PARAMETERS

Subject_Index = {'11', '14', '15'};

%Subject_Index = {'09'};
%% read exp.xlsx DAY1

for si=1:length(Subject_Index)

    % Read xlsx-file (identical log file saved as .xlsx for easier import using readtable() )
    Rlogfile = readtable([datapath filesep 'ps_' Subject_Index{si} '_day1.xlsx']); % read in responses

    Event = Rlogfile(:,3);
    Codes = Rlogfile(:, 4);

    Event = table2cell(Event); % extract event types
    Codes = table2array(Codes); % extract codes aka. response values

    % Find ratings and create list
    start = find(contains (Codes, 'exp_fam2'));
    range_file_C = Codes(start:end, :);
    range_file_E = Event(start:end);
    response = find(contains (range_file_E, 'Response')); % resp is an index (in the end a list of sig. rows) used to search all rows and find all elements in Event which contain the word 'Response'
    resp = range_file_C(response,:); % take the list which indexes rows containing 'Response' to extract all relevant rows from Codes
    resp = resp(3:end-1, :); % delete the first 10 ratings
                              % (4x 1 = enter to continue, 4x training trials, 1x 1 = enter from exp_fam2, 1x repetition of rating from last training trial)
                              % delete the very last rating which is a
                              % repetition of the real last rating

    %% read stimlist.log DAY1

    % Read logfile
    STlogfile = readtable([stimpath filesep 'ps_' Subject_Index{si} '_day1.txt' ]); % read in trials


    stim = STlogfile(10:end-1, 3:6); % the rest is not relevant
    stim = table2cell(stim);

    %% put together
    stim_resp = [stim, resp];
    stim_resp = str2double(string(stim_resp));

    %% change all 2 to 0 and all 3 to 1
    stim_resp(:,5) = changem(stim_resp(:,5), 0, 2); % B = changem(A,new,old); 2 = 0 and 3 = 1
    stim_resp(:,5) = changem(stim_resp(:,5), 1, 3);


    %% save stimuli_response array into .mat
    save([datapath filesep 'ps_' Subject_Index{si} '_day1_stim_resp.mat']);

end

%% read exp.xlsx DAY2
Subject_Index = { '11', '14', '15'};

for si=1:length(Subject_Index)

    % Read xlsx-file (identical log file saved as .xlsx for easier import using readtable() )
    Rlogfile = readtable([datapath filesep 'ps_' Subject_Index{si} '_day2.xlsx']); % read in responses

    Event = Rlogfile(:,3);
    Codes = Rlogfile(:, 4);

    Event = table2cell(Event); % extract event types
    Codes = table2array(Codes); % extract codes aka. response values

    % Find ratings and create list
     % Find ratings and create list
    start = find(contains (Codes, 'exp_fam2'));
    range_file_C = Codes(start:end, :);
    range_file_E = Event(start:end);
    response = find(contains (range_file_E, 'Response')); % resp is an index (in the end a list of sig. rows) used to search all rows and find all elements in Event which contain the word 'Response'
    resp = range_file_C(response,:); % take the list which indexes rows containing 'Response' to extract all relevant rows from Codes
    resp = resp(3:end-1, :); % delete the first 10 ratings
                              % (4x 1 = enter to continue, 4x training trials, 1x 1 = enter from exp_fam2, 1x repetition of rating from last training trial)
                              % delete the very last rating which is a
                              % repetition of the real last rating

    %% read stimlist.log DAY2

    % Read logfile
    STlogfile = readtable([stimpath filesep 'ps_' Subject_Index{si} '_day2.txt' ]); % read in trials


    stim = STlogfile(10:end-1, 3:6); % the rest is not relevant
    stim = table2cell(stim);

    %% put together
    stim_resp = [stim, resp];
    stim_resp = str2double(string(stim_resp));

    %% change all 2 to 0 and all 3 to 1
    stim_resp(:,5) = changem(stim_resp(:,5), 0, 2); % B = changem(A,new,old); 2 = 0 and 3 = 1
    stim_resp(:,5) = changem(stim_resp(:,5), 1, 3);


    %% save stimuli_response array into .mat
    save([datapath filesep 'ps_' Subject_Index{si} '_day2_stim_resp.mat']);

end

%% MLE starts here


    ST = stim_resp(:,1:4);
    R = stim_resp(:,5);
    ModelDimen = 2;
    n_iterations = 5;

    % run 100 iterations with random start values
    % if current one is better (i.e. lower negative log likelihood), keep
    % it
    % otherwise keep the old best one

    best_neg_ll = inf;
    best_position = nan;
    n_problems = 0;

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

    save([datapath filesep 'ps_' Subject_Index{si} '_best_position_day1.mat'], 'best_position');
    save([datapath filesep 'ps_' Subject_Index{si} '_best_position_day2.mat'], 'best_position');
       %To do: Spaces zwischen Subjects alignen
       %procrustes (MATLAB funktion)
       %Für average über Subjects. Procrustes macht spaces so ähnlich wie
       %möglich, ohne die relative position von Stimuli zu verändern
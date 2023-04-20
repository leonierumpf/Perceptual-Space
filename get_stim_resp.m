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
addpath ('C:\Users\rumpf\Documents\perceptual_space\maximum likelihood difference scaling');

%% DEFINE PARAMETERS

Subject_Index = {'02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15'};
d = 1;
%Subject_Index = {'09'};
%Subject_Index = {'10', '11', '12', '13', '14', '15'};
%% read exp.xlsx DAY1
if d == 1
for si=1:length(Subject_Index)

    % Read xlsx-file (identical log file saved as .xlsx for easier import using readtable() )
    Rlogfile = readtable([datapath filesep 'ps_' Subject_Index{si} '_day1.xlsx']); % read in responses

    Event = Rlogfile(:,3);
    Codes = Rlogfile(:, 4);
    Time = Rlogfile(:, 5);
    TTime = Rlogfile(:, 6);

    Event = table2cell(Event); % extract event types
    Codes = table2array(Codes); % extract codes aka. response values
    Time = table2array(Time);
    TTime = table2array(TTime);

    % Find ratings and create list
    start = find(contains (Codes, 'exp_fam2'));
    range_file_C = Codes(start:end, :);
    range_file_E = Event(start:end);
    range_file_T = Time(start:end);
    range_file_TT = TTime(start:end);
   % response = find(contains (range_file_E, 'Response')); % resp is an index (in the end a list of sig. rows) used to search all rows and find all elements in Event which contain the word 'Response'
   % resp = range_file_C(response,:); % take the list which indexes rows containing 'Response' to extract all relevant rows from Codes
   % resp = resp(3:end-1, :); % delete the first 10 ratings
                              % (4x 1 = enter to continue, 4x training trials, 1x 1 = enter from exp_fam2, 1x repetition of rating from last training trial)
                              % delete the very last rating which is a
                              % repetition of the real last rating

                              
    %%%%%%%%%%%%%% Search for all online rating values for BUTTONpress #2
    %%%%%%%%%%%%%% (=upper)
    a =1;
    rating_values_2_all=[];
    rating_values_2_name=[];
    rating_values_2_time=[];
    rating_values_2_reaction=[];

    for p = 1 : length(range_file_C)       %wiederhole solange wie es Reihen im Cellarray Data gibt,
        is_rating_values_2_all = strcmp(range_file_C,'2'); %%% find the buttonpress #7 - stands for NO
        if is_rating_values_2_all(p) > 0
            rating_values_2_all(a,:)    = p;  %%% make a vector with position in the data_logfile
            rating_values_2_name(a,:)   = 00;
            rating_values_2_time(a,1)   = range_file_T(p,:);
            rating_values_2_reaction(a,1)   = range_file_TT(p,:);
            a = a + 1;
        end;
    end;
    
    %%%%%%%%%%%%% Search for all online rating values for BUTTONpress #3
    %%%%%%%%%%%%% (=lower)
    a =1;
    rating_values_3_all=[];
    rating_values_3_name=[];
    rating_values_3_time=[];
    rating_values_3_reaction=[];
    
    for p = 1 : length(range_file_C)       %wiederhole solange wie es Reihen im Cellarray Data gibt,
        is_rating_values_3_all = strcmp(range_file_C,'3'); %%% find the buttonpress #6 - stands for YES
        if is_rating_values_3_all(p) > 0
            rating_values_3_all(a,:)    = p;  %%% make a vector with position in the data_logfile
            rating_values_3_name(a,:)   = 99;
            rating_values_3_time(a,1)   = range_file_T(p,:);
            rating_values_3_reaction(a,1)   = range_file_TT(p,:);
            a = a + 1;
        end;
    end;

    %%%%%%%%%%%%%% Search for rating names (and their order, bc randomized
    %%%%%%%%%%%%%% in the paradigm!)
    a=1;
    D1=[];
    D1_name=[];
    D1_time=[];
    D1_reaction=[];

     for p = 1 : length(range_file_C)       %wiederhole solange wie es Reihen im Cellarray Data gibt,
        is_D1 = strcmp(range_file_C,'D1'); %%% find the string
        if is_D1(p) > 0
            D1(a,:)  = p;  %%% make a vector with position in the logfile
            D1_name(a,:) = 11;  %%% #55 CODE for CS55
            D1_time(a,1) = range_file_T(p,:);
            D1_reaction(a,1) = range_file_TT(p,:);
            a = a + 1;
        end;
    end;

    a=1;
    D2=[];
    D2_name=[];
    D2_time=[];
    D2_reaction=[];

     for p = 1 : length(range_file_C)       %wiederhole solange wie es Reihen im Cellarray Data gibt,
        is_D2 = strcmp(range_file_C,'D2'); %%% find the string
        if is_D2(p) > 0
            D2(a,:)  = p;  %%% make a vector with position in the logfile
            D2_name(a,:) = 22;  %%% #55 CODE for CS55
            D2_time(a,1) = range_file_T(p,:);
            D2_reaction(a,1) = range_file_TT(p,:);
            a = a + 1;
        end;
    end;

    a=1;
    D3=[];
    D3_name=[];
    D3_time=[];
    D3_reaction=[];

     for p = 1 : length(range_file_C)       %wiederhole solange wie es Reihen im Cellarray Data gibt,
        is_D3 = strcmp(range_file_C,'D3'); %%% find the string
        if is_D3(p) > 0
            D3(a,:)  = p;  %%% make a vector with position in the logfile
            D3_name(a,:) = 33;  %%% #55 CODE for CS55
            D3_time(a,1) = range_file_T(p,:);
            D3_reaction(a,1) = range_file_TT(p,:);
            a = a + 1;
        end;
    end;

    a=1;
    D4=[];
    D4_name=[];
    D4_time=[];
    D4_reaction=[];

     for p = 1 : length(range_file_C)       %wiederhole solange wie es Reihen im Cellarray Data gibt,
        is_D4 = strcmp(range_file_C,'D4'); %%% find the string
        if is_D4(p) > 0
            D4(a,:)  = p;  %%% make a vector with position in the logfile
            D4_name(a,:) = 44;  %%% #55 CODE for CS55
            D4_time(a,1) = range_file_T(p,:);
            D4_reaction(a,1) = range_file_TT(p,:);
            a = a + 1;
        end;
    end;

    a=1;
    D5=[];
    D5_name=[];
    D5_time=[];
    D5_reaction=[];

     for p = 1 : length(range_file_C)       %wiederhole solange wie es Reihen im Cellarray Data gibt,
        is_D5 = strcmp(range_file_C,'D5'); %%% find the string
        if is_D5(p) > 0
            D5(a,:)  = p;  %%% make a vector with position in the logfile
            D5_name(a,:) = 55;  %%% #55 CODE for CS55
            D5_time(a,1) = range_file_T(p,:);
            D5_reaction(a,1) = range_file_TT(p,:);
            a = a + 1;
        end;
    end;

    a=1;
    D6=[];
    D6_name=[];
    D6_time=[];
    D6_reaction=[];

     for p = 1 : length(range_file_C)       %wiederhole solange wie es Reihen im Cellarray Data gibt,
        is_D6 = strcmp(range_file_C,'D6'); %%% find the string
        if is_D6(p) > 0
            D6(a,:)  = p;  %%% make a vector with position in the logfile
            D6_name(a,:) = 66;  %%% #55 CODE for CS55
            D6_time(a,1) = range_file_T(p,:);
            D6_reaction(a,1) = range_file_TT(p,:);
            a = a + 1;
        end;
    end;

    a=1;
    D7=[];
    D7_name=[];
    D7_time=[];
    D7_reaction=[];

     for p = 1 : length(range_file_C)       %wiederhole solange wie es Reihen im Cellarray Data gibt,
        is_D7 = strcmp(range_file_C,'D7'); %%% find the string
        if is_D7(p) > 0
            D7(a,:)  = p;  %%% make a vector with position in the logfile
            D7_name(a,:) = 77;  %%% #55 CODE for CS55
            D7_time(a,1) = range_file_T(p,:);
            D7_reaction(a,1) = range_file_TT(p,:);
            a = a + 1;
        end;
    end;

    a=1;
    D8=[];
    D8_name=[];
    D8_time=[];
    D8_reaction=[];

     for p = 1 : length(range_file_C)       %wiederhole solange wie es Reihen im Cellarray Data gibt,
        is_D8 = strcmp(range_file_C,'D8'); %%% find the string
        if is_D8(p) > 0
            D8(a,:)  = p;  %%% make a vector with position in the logfile
            D8_name(a,:) = 88;  %%% #55 CODE for CS55
            D8_time(a,1) = range_file_T(p,:);
            D8_reaction(a,1) = range_file_TT(p,:);
            a = a + 1;
        end;
    end;
    % put all into one array (6name, 6time, 6reaction
%                         7name, 7time, 7 reaction
%                         55name, 55time, 55reaction
%                         145name, 145time, 145reaction)
rating_values_all_time_ordered = [];
rating_values_all_time_ordered(:,1) = [rating_values_2_name;rating_values_3_name;D1_name;D2_name;D3_name;D4_name;D5_name;D6_name;D7_name;D8_name];
rating_values_all_time_ordered(:,2) = [rating_values_2_time;rating_values_3_time;D1_time;D2_time;D3_time;D4_time;D5_time;D6_time;D7_time;D8_time];
rating_values_all_time_ordered(:,3) = [rating_values_2_reaction;rating_values_3_reaction;D1_reaction;D2_reaction;D3_reaction;D4_reaction;D5_reaction;D6_reaction;D7_reaction;D8_reaction];

% sort everything according to time to get right order of stimuli
% presentation
rating_values_all_time_ordered= sortrows(rating_values_all_time_ordered,2);

%%% collect the first button press for CSP- and CSM-stimuli and depending on button make it a yes(6/88=1), no(7/99=0) or missing button (=666)
a=1;b=1;c=1;d=1;e=1;f=1;g=1;h=1;
for ii= 1:length(rating_values_all_time_ordered)
    %%% code 11
    if rating_values_all_time_ordered(ii,1)== 11 && rating_values_all_time_ordered(ii+1,1)== 11 && rating_values_all_time_ordered(ii+2,1)== 00   %upper
        rating(a,1)= 0; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 11 && rating_values_all_time_ordered(ii+1,1)== 11 && rating_values_all_time_ordered(ii+2,1)== 99   %lower
        rating(a,1)= 1; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 11 && rating_values_all_time_ordered(ii+1,1)== 11 && rating_values_all_time_ordered(ii+2,1)~= 00   % no press
        rating(a,1)= 666; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 11 && rating_values_all_time_ordered(ii+1,1)== 11 && rating_values_all_time_ordered(ii+2,1)~= 99   % no press
        rating(a,1)= 666; a=a+1;


    %%% code 22
    elseif rating_values_all_time_ordered(ii,1)== 22 && rating_values_all_time_ordered(ii+1,1)== 22 && rating_values_all_time_ordered(ii+2,1)== 00   %upper
        rating(a,1)= 0; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 22 && rating_values_all_time_ordered(ii+1,1)== 22 && rating_values_all_time_ordered(ii+2,1)== 99   %lower
        rating(a,1)= 1; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 22 && rating_values_all_time_ordered(ii+1,1)== 22 && rating_values_all_time_ordered(ii+2,1)~= 00   % no press
        rating(a,1)= 666; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 22 && rating_values_all_time_ordered(ii+1,1)== 22 && rating_values_all_time_ordered(ii+2,1)~= 99   % no press
        rating(a,1)= 666; a=a+1;

    %%% code 33
   elseif rating_values_all_time_ordered(ii,1)== 33 && rating_values_all_time_ordered(ii+1,1)== 33 && rating_values_all_time_ordered(ii+2,1)== 00   %upper
        rating(a,1)= 0; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 33 && rating_values_all_time_ordered(ii+1,1)== 33 && rating_values_all_time_ordered(ii+2,1)== 99   %lower
        rating(a,1)= 1; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 33 && rating_values_all_time_ordered(ii+1,1)== 33 && rating_values_all_time_ordered(ii+2,1)~= 00   % no press
        rating(a,1)= 666; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 33 && rating_values_all_time_ordered(ii+1,1)== 33 && rating_values_all_time_ordered(ii+2,1)~= 99   % no press
        rating(a,1)= 666; a=a+1;

   %%% code 44
   elseif rating_values_all_time_ordered(ii,1)== 44 && rating_values_all_time_ordered(ii+1,1)== 44 && rating_values_all_time_ordered(ii+2,1)== 00   %upper
        rating(a,1)= 0; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 44 && rating_values_all_time_ordered(ii+1,1)== 44 && rating_values_all_time_ordered(ii+2,1)== 99   %lower
        rating(a,1)= 1; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 44 && rating_values_all_time_ordered(ii+1,1)== 44 && rating_values_all_time_ordered(ii+2,1)~= 00   % no press
        rating(a,1)= 666; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 44 && rating_values_all_time_ordered(ii+1,1)== 44 && rating_values_all_time_ordered(ii+2,1)~= 99   % no press
        rating(a,1)= 666; a=a+1;

   %%% code 55
   elseif rating_values_all_time_ordered(ii,1)== 55 && rating_values_all_time_ordered(ii+1,1)== 55 && rating_values_all_time_ordered(ii+2,1)== 00   %upper
        rating(a,1)= 0; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 55 && rating_values_all_time_ordered(ii+1,1)== 55 && rating_values_all_time_ordered(ii+2,1)== 99   %lower
        rating(a,1)= 1; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 55 && rating_values_all_time_ordered(ii+1,1)== 55 && rating_values_all_time_ordered(ii+2,1)~= 00   % no press
        rating(a,1)= 666; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 55 && rating_values_all_time_ordered(ii+1,1)== 55 && rating_values_all_time_ordered(ii+2,1)~= 99   % no press
        rating(a,1)= 666; a=a+1;

         %%% code 66
   elseif rating_values_all_time_ordered(ii,1)== 66 && rating_values_all_time_ordered(ii+1,1)== 66 && rating_values_all_time_ordered(ii+2,1)== 00   %upper
        rating(a,1)= 0; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 66 && rating_values_all_time_ordered(ii+1,1)== 66 && rating_values_all_time_ordered(ii+2,1)== 99   %lower
        rating(a,1)= 1; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 66 && rating_values_all_time_ordered(ii+1,1)== 66 && rating_values_all_time_ordered(ii+2,1)~= 00   % no press
        rating(a,1)= 666; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 66 && rating_values_all_time_ordered(ii+1,1)== 66 && rating_values_all_time_ordered(ii+2,1)~= 99   % no press
        rating(a,1)= 666; a=a+1;

         %%% code 77
   elseif rating_values_all_time_ordered(ii,1)== 77 && rating_values_all_time_ordered(ii+1,1)== 77 && rating_values_all_time_ordered(ii+2,1)== 00   %upper
        rating(a,1)= 0; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 77 && rating_values_all_time_ordered(ii+1,1)== 77 && rating_values_all_time_ordered(ii+2,1)== 99   %lower
        rating(a,1)= 1; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 77 && rating_values_all_time_ordered(ii+1,1)== 77 && rating_values_all_time_ordered(ii+2,1)~= 00   % no press
        rating(a,1)= 666; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 77 && rating_values_all_time_ordered(ii+1,1)== 77 && rating_values_all_time_ordered(ii+2,1)~= 99   % no press
        rating(a,1)= 666; a=a+1;

         %%% code 88
   elseif rating_values_all_time_ordered(ii,1)== 88 && rating_values_all_time_ordered(ii+1,1)== 88 && rating_values_all_time_ordered(ii+2,1)== 00   %upper
        rating(a,1)= 0; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 88 && rating_values_all_time_ordered(ii+1,1)== 88 && rating_values_all_time_ordered(ii+2,1)== 99   %lower
        rating(a,1)= 1; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 88 && rating_values_all_time_ordered(ii+1,1)== 88 && rating_values_all_time_ordered(ii+2,1)~= 00   % no press
        rating(a,1)= 666; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 88 && rating_values_all_time_ordered(ii+1,1)== 88 && rating_values_all_time_ordered(ii+2,1)~= 99   % no press
        rating(a,1)= 666; a=a+1;

    end
end
    %% read stimlist.log DAY1

    % Read logfile
    STlogfile = readtable([stimpath filesep 'ps_' Subject_Index{si} '_day1.txt' ]); % read in trials


    stim = STlogfile(10:end-1, 3:6); % the rest is not relevant
    stim = table2array(stim);

    %% put together
    stim_resp = [stim, rating];
    stim_resp = str2double(string(stim_resp));
    
    z = 1;
    for idx = 1:length(stim_resp)
        if stim_resp(idx, 5) == 666
            stim_resp(idx,:) = [];
            fprintf('Deleted row %d because of missing value.', idx);
            idx = 1:length(stim_resp);
            z = z + 1;
        end
    end
    %% change all 2 to 0 and all 3 to 1
    %stim_resp(:,5) = changem(stim_resp(:,5), 0, 2); % B = changem(A,new,old); 2 = 0 and 3 = 1
    %stim_resp(:,5) = changem(stim_resp(:,5), 1, 3);


    %% save stimuli_response array into .mat
    save([datapath filesep 'ps_' Subject_Index{si} '_day1_stim_resp.mat']);

end
end
%% read exp.xlsx DAY2

Subject_Index = {'02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15'};
Subject_Index = {'07', '08', '09', '10', '11', '12', '13', '14', '15'};
if d == 2 
for si=1:length(Subject_Index)

     % Read xlsx-file (identical log file saved as .xlsx for easier import using readtable() )
    Rlogfile = readtable([datapath filesep 'ps_' Subject_Index{si} '_day2.xlsx']); % read in responses

    Event = Rlogfile(:,3);
    Codes = Rlogfile(:, 4);
    Time = Rlogfile(:, 5);
    TTime = Rlogfile(:, 6);

    Event = table2cell(Event); % extract event types
    Codes = table2array(Codes); % extract codes aka. response values
    Time = table2array(Time);
    TTime = table2array(TTime);

    % Find ratings and create list
    start = find(contains (Codes, 'exp_fam2'));
    range_file_C = Codes(start:end, :);
    range_file_E = Event(start:end);
    range_file_T = Time(start:end);
    range_file_TT = TTime(start:end);
   % response = find(contains (range_file_E, 'Response')); % resp is an index (in the end a list of sig. rows) used to search all rows and find all elements in Event which contain the word 'Response'
   % resp = range_file_C(response,:); % take the list which indexes rows containing 'Response' to extract all relevant rows from Codes
   % resp = resp(3:end-1, :); % delete the first 10 ratings
                              % (4x 1 = enter to continue, 4x training trials, 1x 1 = enter from exp_fam2, 1x repetition of rating from last training trial)
                              % delete the very last rating which is a
                              % repetition of the real last rating

                              
    %%%%%%%%%%%%%% Search for all online rating values for BUTTONpress #2
    %%%%%%%%%%%%%% (=upper)
    a =1;
    rating_values_2_all=[];
    rating_values_2_name=[];
    rating_values_2_time=[];
    rating_values_2_reaction=[];

    for p = 1 : length(range_file_C)       %wiederhole solange wie es Reihen im Cellarray Data gibt,
        is_rating_values_2_all = strcmp(range_file_C,'2'); %%% find the buttonpress #7 - stands for NO
        if is_rating_values_2_all(p) > 0
            rating_values_2_all(a,:)    = p;  %%% make a vector with position in the data_logfile
            rating_values_2_name(a,:)   = 00;
            rating_values_2_time(a,1)   = range_file_T(p,:);
            rating_values_2_reaction(a,1)   = range_file_TT(p,:);
            a = a + 1;
        end;
    end;
    
    %%%%%%%%%%%%% Search for all online rating values for BUTTONpress #3
    %%%%%%%%%%%%% (=lower)
    a =1;
    rating_values_3_all=[];
    rating_values_3_name=[];
    rating_values_3_time=[];
    rating_values_3_reaction=[];
    
    for p = 1 : length(range_file_C)       %wiederhole solange wie es Reihen im Cellarray Data gibt,
        is_rating_values_3_all = strcmp(range_file_C,'3'); %%% find the buttonpress #6 - stands for YES
        if is_rating_values_3_all(p) > 0
            rating_values_3_all(a,:)    = p;  %%% make a vector with position in the data_logfile
            rating_values_3_name(a,:)   = 99;
            rating_values_3_time(a,1)   = range_file_T(p,:);
            rating_values_3_reaction(a,1)   = range_file_TT(p,:);
            a = a + 1;
        end;
    end;

    %%%%%%%%%%%%%% Search for rating names (and their order, bc randomized
    %%%%%%%%%%%%%% in the paradigm!)
    a=1;
    D1=[];
    D1_name=[];
    D1_time=[];
    D1_reaction=[];

     for p = 1 : length(range_file_C)       %wiederhole solange wie es Reihen im Cellarray Data gibt,
        is_D1 = strcmp(range_file_C,'D1'); %%% find the string
        if is_D1(p) > 0
            D1(a,:)  = p;  %%% make a vector with position in the logfile
            D1_name(a,:) = 11;  %%% #55 CODE for CS55
            D1_time(a,1) = range_file_T(p,:);
            D1_reaction(a,1) = range_file_TT(p,:);
            a = a + 1;
        end;
    end;

    a=1;
    D2=[];
    D2_name=[];
    D2_time=[];
    D2_reaction=[];

     for p = 1 : length(range_file_C)       %wiederhole solange wie es Reihen im Cellarray Data gibt,
        is_D2 = strcmp(range_file_C,'D2'); %%% find the string
        if is_D2(p) > 0
            D2(a,:)  = p;  %%% make a vector with position in the logfile
            D2_name(a,:) = 22;  %%% #55 CODE for CS55
            D2_time(a,1) = range_file_T(p,:);
            D2_reaction(a,1) = range_file_TT(p,:);
            a = a + 1;
        end;
    end;

    a=1;
    D3=[];
    D3_name=[];
    D3_time=[];
    D3_reaction=[];

     for p = 1 : length(range_file_C)       %wiederhole solange wie es Reihen im Cellarray Data gibt,
        is_D3 = strcmp(range_file_C,'D3'); %%% find the string
        if is_D3(p) > 0
            D3(a,:)  = p;  %%% make a vector with position in the logfile
            D3_name(a,:) = 33;  %%% #55 CODE for CS55
            D3_time(a,1) = range_file_T(p,:);
            D3_reaction(a,1) = range_file_TT(p,:);
            a = a + 1;
        end;
    end;

    a=1;
    D4=[];
    D4_name=[];
    D4_time=[];
    D4_reaction=[];

     for p = 1 : length(range_file_C)       %wiederhole solange wie es Reihen im Cellarray Data gibt,
        is_D4 = strcmp(range_file_C,'D4'); %%% find the string
        if is_D4(p) > 0
            D4(a,:)  = p;  %%% make a vector with position in the logfile
            D4_name(a,:) = 44;  %%% #55 CODE for CS55
            D4_time(a,1) = range_file_T(p,:);
            D4_reaction(a,1) = range_file_TT(p,:);
            a = a + 1;
        end;
    end;

    a=1;
    D5=[];
    D5_name=[];
    D5_time=[];
    D5_reaction=[];

     for p = 1 : length(range_file_C)       %wiederhole solange wie es Reihen im Cellarray Data gibt,
        is_D5 = strcmp(range_file_C,'D5'); %%% find the string
        if is_D5(p) > 0
            D5(a,:)  = p;  %%% make a vector with position in the logfile
            D5_name(a,:) = 55;  %%% #55 CODE for CS55
            D5_time(a,1) = range_file_T(p,:);
            D5_reaction(a,1) = range_file_TT(p,:);
            a = a + 1;
        end;
    end;

    a=1;
    D6=[];
    D6_name=[];
    D6_time=[];
    D6_reaction=[];

     for p = 1 : length(range_file_C)       %wiederhole solange wie es Reihen im Cellarray Data gibt,
        is_D6 = strcmp(range_file_C,'D6'); %%% find the string
        if is_D6(p) > 0
            D6(a,:)  = p;  %%% make a vector with position in the logfile
            D6_name(a,:) = 66;  %%% #55 CODE for CS55
            D6_time(a,1) = range_file_T(p,:);
            D6_reaction(a,1) = range_file_TT(p,:);
            a = a + 1;
        end;
    end;

    a=1;
    D7=[];
    D7_name=[];
    D7_time=[];
    D7_reaction=[];

     for p = 1 : length(range_file_C)       %wiederhole solange wie es Reihen im Cellarray Data gibt,
        is_D7 = strcmp(range_file_C,'D7'); %%% find the string
        if is_D7(p) > 0
            D7(a,:)  = p;  %%% make a vector with position in the logfile
            D7_name(a,:) = 77;  %%% #55 CODE for CS55
            D7_time(a,1) = range_file_T(p,:);
            D7_reaction(a,1) = range_file_TT(p,:);
            a = a + 1;
        end;
    end;

    a=1;
    D8=[];
    D8_name=[];
    D8_time=[];
    D8_reaction=[];

     for p = 1 : length(range_file_C)       %wiederhole solange wie es Reihen im Cellarray Data gibt,
        is_D8 = strcmp(range_file_C,'D8'); %%% find the string
        if is_D8(p) > 0
            D8(a,:)  = p;  %%% make a vector with position in the logfile
            D8_name(a,:) = 88;  %%% #55 CODE for CS55
            D8_time(a,1) = range_file_T(p,:);
            D8_reaction(a,1) = range_file_TT(p,:);
            a = a + 1;
        end;
    end;
    % put all into one array (6name, 6time, 6reaction
%                         7name, 7time, 7 reaction
%                         55name, 55time, 55reaction
%                         145name, 145time, 145reaction)
rating_values_all_time_ordered = [];
rating_values_all_time_ordered(:,1) = [rating_values_2_name;rating_values_3_name;D1_name;D2_name;D3_name;D4_name;D5_name;D6_name;D7_name;D8_name];
rating_values_all_time_ordered(:,2) = [rating_values_2_time;rating_values_3_time;D1_time;D2_time;D3_time;D4_time;D5_time;D6_time;D7_time;D8_time];
rating_values_all_time_ordered(:,3) = [rating_values_2_reaction;rating_values_3_reaction;D1_reaction;D2_reaction;D3_reaction;D4_reaction;D5_reaction;D6_reaction;D7_reaction;D8_reaction];

% sort everything according to time to get right order of stimuli
% presentation
rating_values_all_time_ordered= sortrows(rating_values_all_time_ordered,2);

%%% collect the first button press for CSP- and CSM-stimuli and depending on button make it a yes(6/88=1), no(7/99=0) or missing button (=666)
a=1;b=1;c=1;d=1;e=1;f=1;g=1;h=1;
for ii= 1:length(rating_values_all_time_ordered)
    %%% code 11
    if rating_values_all_time_ordered(ii,1)== 11 && rating_values_all_time_ordered(ii+1,1)== 11 && rating_values_all_time_ordered(ii+2,1)== 00   %upper
        rating(a,1)= 0; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 11 && rating_values_all_time_ordered(ii+1,1)== 11 && rating_values_all_time_ordered(ii+2,1)== 99   %lower
        rating(a,1)= 1; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 11 && rating_values_all_time_ordered(ii+1,1)== 11 && rating_values_all_time_ordered(ii+2,1)~= 00   % no press
        rating(a,1)= 666; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 11 && rating_values_all_time_ordered(ii+1,1)== 11 && rating_values_all_time_ordered(ii+2,1)~= 99   % no press
        rating(a,1)= 666; a=a+1;


    %%% code 22
    elseif rating_values_all_time_ordered(ii,1)== 22 && rating_values_all_time_ordered(ii+1,1)== 22 && rating_values_all_time_ordered(ii+2,1)== 00   %upper
        rating(a,1)= 0; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 22 && rating_values_all_time_ordered(ii+1,1)== 22 && rating_values_all_time_ordered(ii+2,1)== 99   %lower
        rating(a,1)= 1; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 22 && rating_values_all_time_ordered(ii+1,1)== 22 && rating_values_all_time_ordered(ii+2,1)~= 00   % no press
        rating(a,1)= 666; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 22 && rating_values_all_time_ordered(ii+1,1)== 22 && rating_values_all_time_ordered(ii+2,1)~= 99   % no press
        rating(a,1)= 666; a=a+1;

    %%% code 33
   elseif rating_values_all_time_ordered(ii,1)== 33 && rating_values_all_time_ordered(ii+1,1)== 33 && rating_values_all_time_ordered(ii+2,1)== 00   %upper
        rating(a,1)= 0; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 33 && rating_values_all_time_ordered(ii+1,1)== 33 && rating_values_all_time_ordered(ii+2,1)== 99   %lower
        rating(a,1)= 1; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 33 && rating_values_all_time_ordered(ii+1,1)== 33 && rating_values_all_time_ordered(ii+2,1)~= 00   % no press
        rating(a,1)= 666; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 33 && rating_values_all_time_ordered(ii+1,1)== 33 && rating_values_all_time_ordered(ii+2,1)~= 99   % no press
        rating(a,1)= 666; a=a+1;

   %%% code 44
   elseif rating_values_all_time_ordered(ii,1)== 44 && rating_values_all_time_ordered(ii+1,1)== 44 && rating_values_all_time_ordered(ii+2,1)== 00   %upper
        rating(a,1)= 0; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 44 && rating_values_all_time_ordered(ii+1,1)== 44 && rating_values_all_time_ordered(ii+2,1)== 99   %lower
        rating(a,1)= 1; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 44 && rating_values_all_time_ordered(ii+1,1)== 44 && rating_values_all_time_ordered(ii+2,1)~= 00   % no press
        rating(a,1)= 666; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 44 && rating_values_all_time_ordered(ii+1,1)== 44 && rating_values_all_time_ordered(ii+2,1)~= 99   % no press
        rating(a,1)= 666; a=a+1;

   %%% code 55
   elseif rating_values_all_time_ordered(ii,1)== 55 && rating_values_all_time_ordered(ii+1,1)== 55 && rating_values_all_time_ordered(ii+2,1)== 00   %upper
        rating(a,1)= 0; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 55 && rating_values_all_time_ordered(ii+1,1)== 55 && rating_values_all_time_ordered(ii+2,1)== 99   %lower
        rating(a,1)= 1; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 55 && rating_values_all_time_ordered(ii+1,1)== 55 && rating_values_all_time_ordered(ii+2,1)~= 00   % no press
        rating(a,1)= 666; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 55 && rating_values_all_time_ordered(ii+1,1)== 55 && rating_values_all_time_ordered(ii+2,1)~= 99   % no press
        rating(a,1)= 666; a=a+1;

         %%% code 66
   elseif rating_values_all_time_ordered(ii,1)== 66 && rating_values_all_time_ordered(ii+1,1)== 66 && rating_values_all_time_ordered(ii+2,1)== 00   %upper
        rating(a,1)= 0; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 66 && rating_values_all_time_ordered(ii+1,1)== 66 && rating_values_all_time_ordered(ii+2,1)== 99   %lower
        rating(a,1)= 1; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 66 && rating_values_all_time_ordered(ii+1,1)== 66 && rating_values_all_time_ordered(ii+2,1)~= 00   % no press
        rating(a,1)= 666; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 66 && rating_values_all_time_ordered(ii+1,1)== 66 && rating_values_all_time_ordered(ii+2,1)~= 99   % no press
        rating(a,1)= 666; a=a+1;

         %%% code 77
   elseif rating_values_all_time_ordered(ii,1)== 77 && rating_values_all_time_ordered(ii+1,1)== 77 && rating_values_all_time_ordered(ii+2,1)== 00   %upper
        rating(a,1)= 0; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 77 && rating_values_all_time_ordered(ii+1,1)== 77 && rating_values_all_time_ordered(ii+2,1)== 99   %lower
        rating(a,1)= 1; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 77 && rating_values_all_time_ordered(ii+1,1)== 77 && rating_values_all_time_ordered(ii+2,1)~= 00   % no press
        rating(a,1)= 666; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 77 && rating_values_all_time_ordered(ii+1,1)== 77 && rating_values_all_time_ordered(ii+2,1)~= 99   % no press
        rating(a,1)= 666; a=a+1;

         %%% code 88
   elseif rating_values_all_time_ordered(ii,1)== 88 && rating_values_all_time_ordered(ii+1,1)== 88 && rating_values_all_time_ordered(ii+2,1)== 00   %upper
        rating(a,1)= 0; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 88 && rating_values_all_time_ordered(ii+1,1)== 88 && rating_values_all_time_ordered(ii+2,1)== 99   %lower
        rating(a,1)= 1; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 88 && rating_values_all_time_ordered(ii+1,1)== 88 && rating_values_all_time_ordered(ii+2,1)~= 00   % no press
        rating(a,1)= 666; a=a+1;
    elseif rating_values_all_time_ordered(ii,1)== 88 && rating_values_all_time_ordered(ii+1,1)== 88 && rating_values_all_time_ordered(ii+2,1)~= 99   % no press
        rating(a,1)= 666; a=a+1;

    end
end

    %% read stimlist.log DAY2

    % Read logfile
    STlogfile = readtable([stimpath filesep 'ps_' Subject_Index{si} '_day2.txt' ]); % read in trials


    stim = STlogfile(10:end-1, 3:6); % the rest is not relevant
    stim = table2array(stim);

    %% put together
    stim_resp = [stim, rating];
    stim_resp = str2double(string(stim_resp));

    z = 1;
    for idx = 1:length(stim_resp)
        if stim_resp(idx, 5) == 666
            stim_resp(idx,:) = [];
            fprintf('Deleted row %d because of missing value (sub%d).', idx, Subject_Index);
            z = z + 1;
        end
    end
    %% change all 2 to 0 and all 3 to 1
    %stim_resp(:,5) = changem(stim_resp(:,5), 0, 2); % B = changem(A,new,old); 2 = 0 and 3 = 1
    %stim_resp(:,5) = changem(stim_resp(:,5), 1, 3);


    %% save stimuli_response array into .mat
    save([datapath filesep 'ps_' Subject_Index{si} '_day2_stim_resp.mat']);

end
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
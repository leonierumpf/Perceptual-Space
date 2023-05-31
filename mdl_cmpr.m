%%% MODEL COMPARISON (AIC)
clear all
close all

clc


%% SET PATH
%
mainpath = 'C:\Users\rumpf\Documents\perceptual_space';
datapath = 'C:\Users\rumpf\Documents\perceptual_space\best_positions';

% Set working dictonary
cd(mainpath);
% all function scripts are stored here:
addpath ('C:\Users\rumpf\Documents\perceptual_space\MLDS scripts');

%% DEFINE PARAMETERS
% these are subject to change

Subject_Index = {'02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15'};

dim1D = 1;
dim2D = 0;

alldays = 0;
dayOne = 0;
dayTwo = 1;

% 1D
if dim1D == 1
    if alldays == 1

        for i = 1:length(Subject_Index)
            load([datapath filesep 'best_positions_alldays1D' filesep 'ps_' Subject_Index{i} '_1D_alldays'], 'MLDS');

            negLL = MLDS.negLL;
            n_params = MLDS.n_params;
            aic = get_aic(negLL, n_params);

            MLDS.AIC = aic;

            save ([datapath filesep 'best_positions_alldays1D' filesep 'ps_' Subject_Index{i} '_1D_alldays'], 'MLDS');

        end

    elseif dayOne == 1
        for i = 1:length(Subject_Index)
            load([datapath filesep 'best_positions_dayOne1D' filesep 'ps_' Subject_Index{i} '_1D_dayOne'], 'MLDS');

            negLL = MLDS.negLL;
            n_params = MLDS.n_params;
            aic = get_aic(negLL, n_params);

            MLDS.AIC = aic;

            save ([datapath filesep 'best_positions_dayOne1D' filesep 'ps_' Subject_Index{i} '_1D_dayOne'], 'MLDS');

        end

    elseif dayTwo == 1
        for i = 1:length(Subject_Index)
            load([datapath filesep 'best_positions_dayTwo1D' filesep 'ps_' Subject_Index{i} '_1D_dayTwo'], 'MLDS');

            negLL = MLDS.negLL;
            n_params = MLDS.n_params;
            aic = get_aic(negLL, n_params);

            MLDS.AIC = aic;

            save ([datapath filesep 'best_positions_dayTwo1D' filesep 'ps_' Subject_Index{i} '_1D_dayTwo'], 'MLDS');

        end
    end
end

% 2D
if dim2D == 1
    if alldays == 1

        for i = 1:length(Subject_Index)
            load([datapath filesep 'best_positions_alldays2D' filesep 'ps_' Subject_Index{i} '_2D_alldays'], 'MLDS');

            negLL = MLDS.negLL;
            n_params = MLDS.n_params;
            aic = get_aic(negLL, n_params);

            MLDS.AIC = aic;

            save ([datapath filesep 'best_positions_alldays2D' filesep 'ps_' Subject_Index{i} '_2D_alldays'], 'MLDS');

        end

    elseif dayOne == 1
        for i = 1:length(Subject_Index)
            load([datapath filesep 'best_positions_dayOne2D' filesep 'ps_' Subject_Index{i} '_2D_dayOne'], 'MLDS');

            negLL = MLDS.negLL;
            n_params = MLDS.n_params;
            aic = get_aic(negLL, n_params);

            MLDS.AIC = aic;

            save ([datapath filesep 'best_positions_dayOne2D' filesep 'ps_' Subject_Index{i} '_2D_dayOne'], 'MLDS');

        end

    elseif dayTwo == 1
        for i = 1:length(Subject_Index)
            load([datapath filesep 'best_positions_dayTwo2D' filesep 'ps_' Subject_Index{i} '_2D_dayTwo'], 'MLDS');

            negLL = MLDS.negLL;
            n_params = MLDS.n_params;
            aic = get_aic(negLL, n_params);

            MLDS.AIC = aic;

            save ([datapath filesep 'best_positions_dayTwo2D' filesep 'ps_' Subject_Index{i} '_2D_dayTwo'], 'MLDS');

        end
    end
end

%%
aic_all_1D = [];
aic_One_1D = [];
aic_Two_1D = [];
aic_all_2D = [];
aic_One_2D = [];
aic_Two_2D = [];

for i  = 1:length(Subject_Index)
    load([datapath filesep 'best_positions_alldays1D' filesep 'ps_' Subject_Index{i} '_1D_alldays'], 'MLDS');
   % load([datapath filesep 'best_positions_dayOne1D' filesep 'ps_' Subject_Index{i} '_1D_dayOne'], 'MLDS');
   % load([datapath filesep 'best_positions_dayTwo1D' filesep 'ps_' Subject_Index{i} '_1D_dayTwo'], 'MLDS');
   % load([datapath filesep 'best_positions_alldays2D' filesep 'ps_' Subject_Index{i} '_2D_alldays'], 'MLDS');
   % load([datapath filesep 'best_positions_dayOne2D' filesep 'ps_' Subject_Index{i} '_2D_dayOne'], 'MLDS');
   % load([datapath filesep 'best_positions_dayTwo2D' filesep 'ps_' Subject_Index{i} '_2D_dayTwo'], 'MLDS');


MLDS.AIC = aic_all_1D(i,1);
% aic = aic_One_1D(i,1);
% aic = aic_Two_1D(i,1);
% aic = aic_all_2D(i,1);
% aic = aic_One_2D(i,1);
% aic = aic_Two_2D(i,1);
save ([datapath filesep 'best_positions_alldays1D'], 'AIC_all_1D');
% save ([datapath filesep 'best_positions_dayOne1D'], 'AIC_One_1D');
% save ([datapath filesep 'best_positions_dayTwo1D'], 'AIC_Two_1D');
% save ([datapath filesep 'best_positions_alldays2D'], 'AIC_all_2D');
% save ([datapath filesep 'best_positions_dayOne2D'], 'AIC_One_2D');
% save ([datapath filesep 'best_positions_dayTwo2D'], 'AIC_Two_2D');
end
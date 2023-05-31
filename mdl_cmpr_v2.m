%%% MODEL COMPARISON (AIC)
clear all
close all

clc


%% SET PATH
%
mainpath = 'C:\Users\rumpf\Documents\perceptual_space';
datapath = 'C:\Users\rumpf\Documents\perceptual_space\best_positions\';

% Set working dictonary
cd(mainpath);
% all function scripts are stored here:
addpath ('C:\Users\rumpf\Documents\perceptual_space\MLDS scripts');
%% DEFINE PARAMETERS
% these are subject to change

Subject_Index = {'02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15'};
bestFITAIC_all = [];
subject_results = randn(6, 3, 14);

for i = 1:length(Subject_Index)
            LL = [];
            n_params = [];

            %model1
            load([datapath filesep 'best_positions_alldays1D' filesep 'ps_' Subject_Index{i} '_1D_alldays'], 'MLDS');

            LL(1,1) = MLDS.negLL*(-1);
            n_params(1,1) = MLDS.n_params - 2;

            %model2
            load([datapath filesep 'best_positions_dayOne1D' filesep 'ps_' Subject_Index{i} '_1D_dayOne'], 'MLDS');

            LL(2,1) = MLDS.negLL*(-1);
            n_params(2,1) = MLDS.n_params - 2;

            %model3
            load([datapath filesep 'best_positions_dayTwo1D' filesep 'ps_' Subject_Index{i} '_1D_dayTwo'], 'MLDS');

            LL(3,1) = MLDS.negLL*(-1);
            n_params(3,1) = MLDS.n_params - 2;

            %model4
            load([datapath filesep 'best_positions_alldays2D' filesep 'ps_' Subject_Index{i} '_2D_alldays'], 'MLDS');

            LL(4,1) = MLDS.negLL*(-1);
            n_params(4,1) = MLDS.n_params - 4;

            %model5
            load([datapath filesep 'best_positions_dayOne2D' filesep 'ps_' Subject_Index{i} '_2D_dayOne'], 'MLDS');

            LL(5,1) = MLDS.negLL*(-1);
            n_params(5,1) = MLDS.n_params - 4;

            %model6
            load([datapath filesep 'best_positions_dayTwo2D' filesep 'ps_' Subject_Index{i} '_2D_dayTwo'], 'MLDS');

            LL(6,1) = MLDS.negLL*(-1);
            n_params(6,1) = MLDS.n_params - 4;

            Tbl = table(LL,n_params,RowNames="Model"+string(1:6));
            aic = aicbic(LL,n_params);

            %matrix mit dimensionen: ModelDimensions x Parameters x Subjects 
            subject_results(:, 1, i) = LL;
            subject_results(:, 2, i) = n_params;
            subject_results(:, 3, i) = aic;

            [~,idxmin] = min(aic);
            bestFitAIC = Tbl.Properties.RowNames{idxmin};

            bestFITAIC_all(i,1) = idxmin;
end

save ([datapath], 'bestFITAIC_all');
%% RETROGRADE MODELLING %%
% This code takes as input the fitted mRNA and the logistic parameters of
% the protein of choice. It calculates the turnover rates for synthesis and
% degradation.

% Solve analytically for M(t) by using as input the experimental protein expression P(t)
% The synthesis and degradation rates Sp and Dp are calculated by non-linear regression 
% between the calculated solution and the corresponding experimental data 

% See for reference Matsubayashi, Y, S�nchez-S�nchez, BJ, et al. Developmental Cell (2020).
% Copyright - Stefania Marcotti (Stramer lab) 2020
% For issues please contact us at https://www.stramerlab.com/contact.html

%% Load input and choose parameters %%

% import fitted mRNA
uiwait(msgbox('Import fitted mRNA (.csv file)'))
[file, d] = uigetfile('*.csv');
data = readtable(fullfile(d,file), 'ReadVariableNames', 0);
data = table2array(data);
data(any(isnan(data), 2), :) = [];
time = data(:,1);
mRNA = data(:,2);

% import logistic parameters for sample of choice
uiwait(msgbox('Import logistic parameters for protein (.csv file)'))
[file_log, d_log] = uigetfile('*.csv');
data_log = readtable(fullfile(d_log, file_log), 'ReadVariableNames', 0);
data_log = table2array(data_log);

% set parameters
prompt = {'Time interval [min]', ...
    'Starting time experiment [h]',...
    'Finishing time experiment [h]'};
title_prompt = 'User parameters';
dims = [1 35];
user_answer = inputdlg(prompt,title_prompt,dims);

% define temporal span
idx = find(time >= str2double(user_answer{2,1}) & time <= str2double(user_answer{3,1}));
time_exp = time(idx);	% [h]
mRNA = mRNA(idx);

% set initial values for fitting
beta0 = [0.01 0.01];

%% Find turnover rates %%

Sp_star = zeros(size(data_log,2),1);
Dp_star = zeros(size(data_log,2),1);
Sp_star_confidence_intervals = zeros(size(data_log,2),2);
Dp_star_confidence_intervals = zeros(size(data_log,2),2);

% for each logistic curve
for q = 1:size(data_log,2)
    
    % P_exp (experimental value)
    K = data_log(1,q);
    r = data_log(3,q);
    ti = data_log(2,q);
    
    % perform fitting
    funct_parameters = [K, r, ti];
    [beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(time_exp, mRNA, @(beta, time_exp)retrograde_funct(beta, funct_parameters, time_exp), beta0);
    ci = nlparci(beta,R,'Jacobian',J);

    % save current values in output vectors
    Sp_star(q,1) = beta(1);
    Dp_star(q,1) = beta(2);
    Sp_star_confidence_intervals(q,:) = ci(1,:);
    Dp_star_confidence_intervals(q,:) = ci(2,:);
    
end

%% Save output as .csv files %%

writetable(array2table(Sp_star), fullfile(d,'retrograde_Sp_star.csv'), 'WriteVariableNames', 0)
writetable(array2table(Dp_star), fullfile(d,'retrograde_Dp_star.csv'), 'WriteVariableNames', 0)
writetable(array2table(Sp_star_confidence_intervals), fullfile(d,'retrograde_Sp_star_confidence_intervals.csv'), 'WriteVariableNames', 0)
writetable(array2table(Dp_star_confidence_intervals), fullfile(d,'retrograde_Dp_star_confidence_intervals.csv'), 'WriteVariableNames', 0)

clear; clc
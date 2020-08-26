%% ANTEROGRADE MODELLING %%
% This code takes as input the fitted mRNA and the logistic parameters of
% the protein of choice. It calculates the turnover rates for synthesis and
% degradation. 
% See for reference Matsubayashi, Y, Sánchez-Sánchez, BJ, et al. Developmental Cell (2020).
% Copyright - Stefania Marcotti (Stramer lab) 2020
% For issues please contact us at https://www.stramerlab.com/contact.html

%% Load input and choose parameters %%

% import fitted mRNA
uiwait(msgbox('Import fitted mRNA (.csv file)'))
[file, d] = uigetfile('*.csv');
data = readtable(fullfile(d,file));
data = table2array(data);
time = data(:,1);
mRNA = data(:,2);

% import logistic parameters for sample of choice
uiwait(msgbox('Import logistic parameters for protein (.csv file)'))
[file_log, d_log] = uigetfile('*.csv');
data_log = readtable(fullfile(d_log, file_log));
data_log = table2array(data_log);

% set parameters
prompt = {'Time interval [min]', ...
    'Starting time experiment [h]',...
    'Finish time experiment [h]'};
title_prompt = 'User parameters';
dims = [1 35];
definput = {'2','-10','30'};
user_answer = inputdlg(prompt,title_prompt,dims,definput);

% define temporal span
delta_t = str2double(user_answer{1,1})/60;	% [h]
time_exp = (str2double(user_answer{2,1}):delta_t:str2double(user_answer{3,1}))';	% [h]
idx = find(time >= min(time_exp) & time <= max(time_exp));
mRNA = mRNA(idx);

% set initial values for fitting
beta0 = [0.01 0.01];

%% Find turnover rates %%

% initialise output
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
    P_exp = K ./ (1+exp(-r .* (time_exp - ti)));
    
    % perform fitting 
    funct_parameters = [delta_t; mRNA];
    [beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(time_exp, P_exp, @(beta, time_exp)anterograde_funct(beta, funct_parameters, time_exp), beta0);
    ci = nlparci(beta,R,'Jacobian',J);
    
    % save current values in output vectors
    Sp_star(q,1) = beta(1);
    Dp_star(q,1) = beta(2);
    Sp_star_confidence_intervals(q,:) = ci(1,:);
    Dp_star_confidence_intervals(q,:) = ci(2,:);

    % clear reused variables
    clear P_exp

end

%% Save output as .csv files %%

writetable(array2table(Sp_star), fullfile(d,'Sp_star.csv'), 'WriteVariableNames', 0)
writetable(array2table(Dp_star), fullfile(d,'Dp_star.csv'), 'WriteVariableNames', 0)
writetable(array2table(Sp_star_confidence_intervals), fullfile(d,'Sp_star_confidence_intervals.csv'), 'WriteVariableNames', 0)
writetable(array2table(Dp_star_confidence_intervals), fullfile(d,'Dp_star_confidence_intervals.csv'), 'WriteVariableNames', 0)

clear; clc
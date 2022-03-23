% clear all
% close all
clc
format long
set(0, ...
    'DefaultLineLineWidth',1.2, ...
    'DefaultAxesFontSize', 18 ...    
);

% Input info
nameDictPath = '.\aux-files\';
Params = jsondecode(fileread(strcat(nameDictPath, 'absorbance.json')));
Params.data_directory = '.\datalake\'; 

% Add paths to Matlab
loadersDir = '.\loaders\';
utisDir = '.\utils\';
dataDetailsDir =  strcat(Params.data_directory, 'data-details\');
addpath(Params.data_directory, loadersDir, nameDictPath, utisDir);

% Read experiment catalog
[~, ~, expCatalog] = xlsread(strcat(dataDetailsDir, 'exp_catalog.csv'));
expCatalog = cell2table(expCatalog(2:end, :), 'VariableNames', expCatalog(1, :))

% The next 2 lines allow interative mode. 
% prompt = 'Give me the experiment code: ';
% expCode = input(prompt, 's');
expCode = 'p000';

files = dir(fullfile(Params.data_directory, strcat(expCode,'*.txt')));
% prompt = 'Do you want to proceed Y/N [Y]: ';
% answer = input(prompt, 's');

filesNames = {files.name};
disp(filesNames)


%% Here the analysis starts calling a class with all the methods to load, transform and show the data.

EMPTY = SPECTRA.empty;
WATER = SPECTRA.empty;
ABSSPECTRA = SPECTRA.empty;

for i = 1:length(filesNames)
    
    fList = strsplit(filesNames{i}, '');
    % fileID = fList{1};
    fileID = lower(filesNames{i}(5));

    switch fileID
        case 'b'
            EMPTY = [EMPTY, SPECTRA(filesNames{i}, Params.data_directory)];
            EMPTY(end).Plot('intensity');
        case 'r'
            WATER = [WATER, SPECTRA(filesNames{i}, Params.data_directory)];
            WATER(end).Plot('intensity');
        case 's'
            ABSSPECTRA = [ABSSPECTRA, SPECTRA(filesNames{i}, Params.data_directory)];
            ABSSPECTRA(end).Plot('intensity');
    end
end

WATER.SetBackground(EMPTY);
ABSSPECTRA.SetBackground(EMPTY);
ABSSPECTRA.SetReference(WATER);
ABSSPECTRA.GetAbsorbance();
ABSSPECTRA.Plot('absorbedPow');
ABSSPECTRA.Plot('absorbance');
ABSSPECTRA.FitAbsorbance(Params.fitting);

%plotLines = flip(get(gca, 'Children'));
% legend(plotLines(1:5), { ...
%     '5 mg/ml', ...
%     'Fitting', ...
%     'Used region' ...    
% }, 'Location', 'best');

%% End of the script
% clear all
% close all
clc
format long

% Input info
loadersDir = '.\loaders\';
utisDir = '.\utils\';
nameDictPath = '.\aux-files\';

% Include the data folder to Matlab's path
addpath(loadersDir, nameDictPath, utisDir);

% Initial parameters
BeamProfParams = jsondecode(fileread(strcat(nameDictPath, 'beamprofile.json')));
BeamProfParams.data_directory = '.\datalake\'; 

% Get data details and experiemnt cathalog
dataDetailsDir =  strcat(BeamProfParams.data_directory, 'data-details\');

[~, ~, sampleCatalog] = xlsread(strcat(dataDetailsDir, 'samples_catalog.csv'));
sampleCatalog = cell2table(sampleCatalog(2:end, :), 'VariableNames', sampleCatalog(1, :));

[~, ~, expCatalog] = xlsread(strcat(dataDetailsDir, 'exp_catalog.csv'));
expCatalog = cell2table(expCatalog(2:end, :), 'VariableNames', expCatalog(1, :))

% Interactive mode commented...
% prompt = 'Give me the experiment code: ';
% expCode = lower(input(prompt, 's'));
% prompt = 'Give me the file name: ';
% file = lower(input(prompt, 's'));
% BeamProfParams.beam_profile.fileName = feval(@(v) v{1}, expCatalog{ contains(lower(expCatalog.File), strcat(file, '.png')), 1});

expCode = 'p001'; % Esperiment code of the picture used to evaluate the beam profile of the laser. 
expCatalog = expCatalog(contains(lower(expCatalog.ExperimentCode), expCode), :)


%% Analysis starts here. 
BEAMPROFILE = BEAM_PROFILE( ...
    BeamProfParams.beam_profile.fileName, ...
    BeamProfParams.data_directory, ...
    BeamProfParams.beam_profile.width_method, ...
    BeamProfParams.beam_profile.px_size ...
);

BEAMPROFILE.RemoveBackground();
BEAMPROFILE.Trim(2); 
BEAMPROFILE.Normalize();
BEAMPROFILE.Width(BEAMPROFILE.camPxSize);

disp('File:')
BEAMPROFILE

%  plots
BEAMPROFILE.Imshow();
BEAMPROFILE.Surf('image');
BEAMPROFILE.Surf('gaussian');

%% End of script
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
loadersDir = '.\loaders\';
utisDir = '.\utils\';

% Initial parameters
Params = jsondecode(fileread(strcat(nameDictPath, 'absorbance.json')));
Params.data_directory = '.\datalake\'; 
dataDetailsDir =  strcat(Params.data_directory, 'data-details\');

% Add path top the Matlab environment
addpath(Params.data_directory, loadersDir, nameDictPath, utisDir);

% Get data details and experiments cathalog
[~, ~, expCatalog] = xlsread(strcat(dataDetailsDir, 'exp_catalog.csv'));
expCatalog = cell2table(expCatalog(2:end, :), 'VariableNames', expCatalog(1, :))

% Interactive mode commented...
% prompt = 'Give me the file name: ';
% expCode = input(prompt, 's');
expCode = 'p000e';

% files = dir(fullfile(Params.data_directory, strcat(expCode,'*.txt')));
% prompt = 'Do you want to proceed Y/N [Y]: ';
% answer = input(prompt, 's');
file = dir(fullfile(Params.data_directory, strcat(expCode,'.csv')));

filesNames = {file.name};
disp(filesNames)


%% Analsys starts here
EMSPECTRA = SPECTRA(filesNames{1}, Params.data_directory);
EMSPECTRA.laserCurrent = '600 nm';
EMSPECTRA.RemoveBaseLine();
EMSPECTRA.NormalizeMod([790, 810]);

% Show results
figure(20)
% SP.Wavelength(), SP.Intensity(), '-m', ...
% SP1.Wavelength(), SP1.Intensity(), '-b', ... 
plot( ... 
    EMSPECTRA.Wavelength(), EMSPECTRA.Intensity(), '-', 'DisplayName', filesNames{1} ... 
)
xlim([740,860])
% ylim([0.001,1.1])
% title(sampleName)
hold on
% legend( ...
%     EMSPECTRA.laserCurrent
% );
xlabel('Wavelength (nm)');
ylabel('Intensity (a.u.)');

%% End of script
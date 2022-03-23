%%%%%%%%%% This is the initial version of the code to analyse the QY data of the UCNPs. 
%%%%%%%%%% This code is used to generate the luminescence curve as an example of possible measurement with QY system.
%%%%%%%%%% However, this code doesn't have all the needed features to compute QY. The most updated code is written in Python in the ../quantum-yield/ folder.
  
% clear all
% close all
clc
format long
set(0, ...
    'DefaultLineLineWidth',1.2, ...
    'DefaultAxesFontSize', 18 ...    
);


% Input info
loadersDir = '.\loaders\';
utisDir = '.\utils\';
nameDictPath = '.\aux-files\';
% Include the data folder to Matlab's path
addpath(loadersDir, nameDictPath, utisDir);

% Initial parameters from json file
Params = jsondecode(fileread(strcat(nameDictPath, 'lumsettings.json')));
Params.data_directory = '.\datalake\'; 

% Load data details
dataDetailsDir =  strcat(Params.data_directory, 'data-details\');

[~, ~, sampleCatalog] = xlsread(strcat(dataDetailsDir, 'samples_catalog.csv'));
sampleCatalog = cell2table(sampleCatalog(2:end, :), 'VariableNames', sampleCatalog(1, :));

[~, ~, expCatalog] = xlsread(strcat(dataDetailsDir, 'exp_catalog.csv'));
expCatalog = cell2table(expCatalog(2:end, :), 'VariableNames', expCatalog(1, :))

% Interactive mode commented
% prompt = 'Give me the experiment code: ';
% expCode = lower(input(prompt, 's'));
expCode = 'p000';
expCatalog = expCatalog(contains(lower(expCatalog.ExperimentCode), expCode), :)

% prompt = 'Give me spot size: ';
% spotSize = lower(input(prompt, 's'));
spotSize = 's1';

[~, ~, absorbanceData] = xlsread(strcat(Params.data_directory, strcat(expCode, '.xlsx')), 'absorp', 'B2:H5');
absorbanceData = cell2table(absorbanceData(2:end, :));

Params.beam_profile.laser976 = feval(@(v) v{1}, expCatalog{ contains(lower(expCatalog.File), strcat('l976', spotSize, '.png')), 1});
Params.beam_profile.laser785 = feval(@(v) v{1}, expCatalog{ contains(lower(expCatalog.File), strcat('l785', spotSize, '.png')), 1});

expCatalog = join(expCatalog(contains(expCatalog.Measurement, 'luminescence'), :), sampleCatalog, 'Keys', 'Sample');

ucnpFileCatalog = expCatalog(contains(lower(expCatalog.File), strcat('u', spotSize, '.csv')), :);

% Updating parameters regarding the UCNPs
Params.ucnp.file_name = ucnpFileCatalog.File{1};
Params.ucnp.quantum_yield = ucnpFileCatalog.QuantumYield;
Params.ucnp.concentration = ucnpFileCatalog.FinalConcentration_mgml{1};
Params.ucnp.filter_transmitance = Params.system.filter_transmitance976;
Params.ucnp.absorbance = absorbanceData{3, 6};
Params.ucnp.absorbance_std = absorbanceData{3, 7};
Params.ucnp.solvent.absorbance = absorbanceData{1, 6}; 
Params.ucnp.solvent.absorbance_std = absorbanceData{1, 7}; 
Params.ucnp.solvent.name = ucnpFileCatalog.Solvent{1};
Params.ucnp.solvent.refractive_idx = feval(@(t) t.refractiveIndex, sampleCatalog(contains(sampleCatalog.Sample, ucnpFileCatalog.Solvent{1}), :));
Params.ucnp.label = feval(@(s) s{1}, strcat(ucnpFileCatalog.ExperimentCode, '.', ucnpFileCatalog.Sample, '-', ucnpFileCatalog.SpotSize));


dyeFileCatalog = expCatalog(contains(lower(expCatalog.File), strcat('d', spotSize, '.csv')), :);

% Updating parameters regarding the Dye
Params.dye.file_name = dyeFileCatalog.File{1};
Params.dye.quantum_yield = dyeFileCatalog.QuantumYield;
Params.dye.concentration = dyeFileCatalog.FinalConcentration_mgml{1};
Params.dye.filter_transmitance = Params.system.filter_transmitance785;
Params.dye.absorbance = absorbanceData{3, 2};
Params.dye.absorbance_std = absorbanceData{3, 3};
Params.dye.solvent.absorbance = absorbanceData{1, 2};
Params.dye.solvent.absorbance_std = absorbanceData{1, 3};
Params.dye.solvent.name = dyeFileCatalog.Solvent{1};
Params.dye.solvent.refractive_idx = feval(@(t) t.refractiveIndex, sampleCatalog(contains(sampleCatalog.Sample, dyeFileCatalog.Solvent{1}), :));
Params.dye.label = feval(@(s) s{1}, strcat(dyeFileCatalog.ExperimentCode, '.', dyeFileCatalog.Sample, '-', dyeFileCatalog.SpotSize));


%% Analysis starts here. 

% Beam profile of the laser 785nm:
BEAMPROFILE785 = BEAM_PROFILE( ...
    Params.beam_profile.laser785, ...
    Params.data_directory, ...
    Params.beam_profile.width_method, ...
    Params.beam_profile.px_size ...
);
BEAMPROFILE785.RemoveBackground();
BEAMPROFILE785.Trim(1); 
BEAMPROFILE785.Normalize();

% Beam profile of laser 976nm:
BEAMPROFILE976 = BEAM_PROFILE( ...
    Params.beam_profile.laser976, ...
    Params.data_directory, ...
    Params.beam_profile.width_method, ...
    Params.beam_profile.px_size ...
);
BEAMPROFILE976.RemoveBackground();
BEAMPROFILE976.Trim(1);
BEAMPROFILE976.Normalize();

% UCNPs and Dye data from QY system loaded here
UCNP = SWEEP_MODE( Params.ucnp, Params.data_directory );
DYE = SWEEP_MODE( Params.dye, Params.data_directory );
    
% Remove background
UCNP.RemoveBackground();
DYE.RemoveBackground();

% Calculates luminescence Signal
UCNP.LuminSignal();
DYE.LuminSignal();

% Obtains absorbed Power at the centre (mW/cm)
UCNP.AbsorbedPowAtCentre();
DYE.AbsorbedPowAtCentre();

% Calculates experimental Quantum Yield
DYE.ExpQuantumYield();

% Relative Quantum Yield
UCNP.RelQuantumYield(DYE);

% Obtains power density
UCNP.FluenceRate(BEAMPROFILE976.Area());
DYE.FluenceRate(BEAMPROFILE785.Area());

% UCNPs and Reference Dye
ucnpLuminFit = polyval( ...
    polyfit(UCNP.absorbedPowAtCentre, UCNP.luminSignal, 3), ...
    UCNP.absorbedPowAtCentre ...
);

refLuminFit = polyval( ...
    polyfit(DYE.absorbedPowAtCentre, DYE.luminSignal, 1), ...
    DYE.absorbedPowAtCentre ...
);

DYE.Plot('luminescence');
figure(104)
plot(DYE.absorbedPowAtCentre, refLuminFit, 'r-', 'DisplayName', 'Fitting');

UCNP.Plot('luminescence', Params.ucnp.label);
figure(104)
plot(UCNP.absorbedPowAtCentre, ucnpLuminFit, 'r-', 'DisplayName', 'Fitting');

% Plotting the results
UCNP.Plot('expQY', Params.ucnp.label);
DYE.Plot('expQY', Params.dye.label);
UCNP.Plot('relQY', Params.ucnp.label);
DYE.Plot('rawdata', Params.dye.label);
UCNP.Plot('rawdata', Params.ucnp.label);

%% End of QY analysis 
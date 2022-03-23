% clear all
close all
clc
format long

% Input info
loadersDir = '.\loaders\';
utisDir = '.\utils\';
nameDictPath = '.\aux-files\';

% Include the data folder to Matlab's path
addpath(loadersDir, nameDictPath, utisDir);

% Run all the modules first, then run this one
%   module               # filename
%   A_beamprofile        # p001l976s1.png
%   B_quantumyield       # p000US1/p000DS1
%   C_spectra            # p000b/s/r.txt
%   D_emissionspectra    # p000e.txt
%   F_autocorrelation    # p000c.txt

%% Plots settings
set(0, ...
    'DefaultLineLineWidth',1.8, ...
    'defaultTextFontName', 'Halvetica', ...
    'DefaultAxesFontSize', 12 ...    
);
colormat = [ ...
    245,178,99; ...
    255,178,38; ...
    232,134,37; ...
    0,100,0; ...
    102,51,153; ...
    0,139,139 ...
] / 255;
%%

%% Prepare data for plotting and making the figures
BEAMPROFILE.Reset()
BEAMPROFILE.RemoveBackground();
centre = double(BEAMPROFILE.MaxPeak());
beamWidthXY = BEAMPROFILE.Width(1); % in pxs
xLimArr1 = [centre(3) - 1.4 * beamWidthXY(1), centre(3) + 1.4 * beamWidthXY(1)];
yLimArr1 = [centre(2) - 1.25 * beamWidthXY(1), centre(2) + 1.25 * beamWidthXY(1)];
BEAMPROFILE.Trim(2); 

% plots Beam profile
screenSize = get(0, 'ScreenSize');
f4 = figure (4);
set(f4, 'position', [screenSize(3) * 0.2, screenSize(4) * 0.05, screenSize(3) * 0.55, screenSize(4) * 0.85]);
numOfRows = 22;
numOfCols = 23;
gridMatrix = vec2mat([1:numOfRows * numOfCols], numOfCols);

% subplot(6, 2, 5)
%     1:10 11:13 13:23 
%     1:10 11:13 13:14 14:23
% 1:2
% 4:10 
% 11:12
% 13:22
%%

%% Figure 4. 
figure (4);
gM = gridMatrix(16:22, 16:23);
gM = sort(gM(:));
subplot(numOfRows, numOfCols, gM)
img = imshow(BEAMPROFILE.rawData);
xlabel('Px', 'Interpreter', 'Latex')
ylh = ylabel('Px', 'Interpreter', 'Latex');
xlim(xLimArr1);
ylim(yLimArr1);
ylh.Position(1) = 955;


BEAMPROFILE.Normalize();
centre = double(BEAMPROFILE.MaxPeak());
beamWidthXY = BEAMPROFILE.Width(1); % in pxs
xLimArr = [- 1.45 * beamWidthXY(1), 1.45 * beamWidthXY(1)];

% top
gM = gridMatrix(13:14, 16:23);
gM = sort(gM(:));
subplot(numOfRows, numOfCols, gM)
pxRange = length(BEAMPROFILE.data(:, centre(2)));
plot([-floor(pxRange / 2) : floor(pxRange / 2)], BEAMPROFILE.data(:, centre(2)), 'k.')
set(gca, 'ytick', []);
xlim(xLimArr);

% left
gM = gridMatrix(16:22, 13:14);
gM = sort(gM(:));
subplot(numOfRows, numOfCols, gM)
plot(-1 * BEAMPROFILE.data(centre(3), :), [-floor(pxRange / 2) : floor(pxRange / 2)], 'k.')
ax = gca;
ax.YLim = xLimArr; 
ax.YAxisLocation = 'right';
ax.XTick = [];
set(gca, 'xtick', []);
ylim(xLimArr);


% plots Luminescence
figure (4);
gM = gridMatrix(1:10, 14:23);
gM = sort(gM(:));
subplot(numOfRows, numOfCols, gM)
plot(UCNP.absorbedPowAtCentre, UCNP.luminSignal, 'ko', 'MarkerSize', 3, 'DisplayName', 'UCNP', 'color', 'black', 'markerfacecolor', [1, 1, 1]);
hold on
plot(UCNP.absorbedPowAtCentre, ucnpLuminFit, 'r-', 'DisplayName', 'Fitting');
xlim([0, 17])
ylim([0, 4])
% title('Luminescence signal')
xlabel('Absorbed power per unit length (mW/cm)', 'Interpreter', 'latex')
ylabel('Luminescence signal (counts)', 'Interpreter', 'latex')

% Absorbance and scattering
[abArr, model] = ABSSPECTRA.FitAbsorbance(Params.fitting);
from = ABSSPECTRA.FindClosest(ABSSPECTRA.Wavelength(), 500);
to = ABSSPECTRA.FindClosest(ABSSPECTRA.Wavelength(), 900);

% plots absorpsion spectra
figure (4);
gM = gridMatrix(1:10, 1:10);
gM = sort(gM(:));
subplot(numOfRows, numOfCols,gM)
semilogy(ABSSPECTRA.Wavelength(), ABSSPECTRA.Absorbance(), '.', 'color', colormat(4, :), 'markerfacecolor', [1, 1, 1]);
hold on
semilogy(ABSSPECTRA.Wavelength(), model(abArr, ABSSPECTRA.Wavelength()), 'k-');
semilogy(ABSSPECTRA.data(from:to, 1), model(abArr, ABSSPECTRA.data(from:to, 1)), 'r-');
set(gca, 'YTick', [10 ^ -1, 2.5 * 10 ^ -1, 5 * 10 ^ -1, 10 ^ 0 ])
% title(sampleName)
xlabel('Wavelength (nm)', 'Interpreter', 'latex');
ylabel('Absorbance (a.u.)', 'Interpreter', 'latex');
xlim([500 1050])
ylim([0.14, 1.3])


% plots Emission spectra
figure (4);
gM = gridMatrix(13:22, 1:10);
gM = sort(gM(:));
subplot(numOfRows, numOfCols, gM)
plot(EMSPECTRA.Wavelength(), EMSPECTRA.Intensity(), '-', 'MarkerSize', 3, 'color', colormat(6, :), 'markerfacecolor', [1, 1, 1]);
xlim([350,1100])
ylim([0, 1])
xlabel('Wavelength (nm)', 'Interpreter', 'latex');
ylabel('Intensity (a.u.)', 'Interpreter', 'latex');

saveas(f4, './plots/Fig4.png')
%%

%% Figure 3
screenSize = get(0, 'ScreenSize');
f3 = figure(3);
set(f3, 'position', [screenSize(3) * 0.2, screenSize(4) * 0.05, screenSize(3) * 0.5, screenSize(4) * 0.85]);
numOfRows = 2;
numOfCols = 2;

subplot(numOfRows, numOfCols, [1]);
s = surf(BEAMPROFILE.data);
colormap(get(s, 'Parent'), parula)
view(90, -90);
box on;
set(s, 'edgecolor','none');
xlabel('Px', 'Interpreter', 'latex')
ylabel('Px', 'Interpreter', 'latex')
dataShape = size(BEAMPROFILE.data);
ylim([0, dataShape(1)]);
xlim([0, dataShape(2)]);

subplot(numOfRows, numOfCols, [2]);
s = surf(BEAMPROFILE.data);
colormap(get(s, 'Parent'), parula)
box on;
set(s, 'edgecolor','none');
xlabel('Px', 'Interpreter', 'latex')
ylabel('Px', 'Interpreter', 'latex')
zlabel('Intensity (a.u)', 'Interpreter', 'latex')
dataShape = size(BEAMPROFILE.data);
ylim([0, dataShape(1)]);
xlim([0, dataShape(2)]);

% plots correlation
subplot(numOfRows, numOfCols, [3]);
semilogx(x, 'color', colormat(6, :));
ylim([1, 2])
xlabel('Time (${\mu}s$)', 'Interpreter', 'latex')
ylabel('APD Signal ($V$)', 'Interpreter', 'latex')
set(gca, 'XTick', [1e0, 1e2, 1e4, 1e6, 1e8 ])

subplot(numOfRows, numOfCols, [4]);
semilogx(f, P1, 'color', colormat(6, :));
xlabel('Frequency (Hz)', 'Interpreter', 'latex')
ylabel('$\left |P(f) \right |$', 'Interpreter', 'latex')
set(gca, 'XTick', [1e-1, 1e1, 1e3, 1e5 ])

saveas(f3, './plots/Fig3.png')

%%
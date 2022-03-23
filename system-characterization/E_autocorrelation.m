close all
% clear all

filePath = '.\datalake\';

% filename = 'background_1Ms_60s_b.txt';
filename = 'p000c.txt'; % 700mA_1Ms_60s_a.txt

data = importdata(strcat(filePath, filename));
% delimiter = {''};
% formatSpec = '%f%[^\n\r]';
% fileID = fopen(filename,'r');
% dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
% fclose(fileID);
% x = table(dataArray{1:end-1}, 'VariableNames', {'VarName1'});
% clearvars filename delimiter formatSpec fileID dataArray ans;

% x = data(1:10000000);
x = data;
Fs = 1.0e+6; % Hz
T = 1/Fs;
L = length(x);
t = (0: L - 1) * T; % s
Y = fft(x); 
P2 = abs(Y/L);
P1 = P2(1:L/2 + 1);
P1(2:end-1) = 2 * P1(2:end-1);

f = Fs*(0:(L/2)) / L;
semilogx(f, P1);
xlabel('Frequency [Hz]', 'Interpreter', 'latex')
ylabel('P(f)', 'Interpreter', 'latex')

% time = linspace(1:)
%load('autocorr')
%measurement_no=1;
% [c,lags]=xcorr(x,'unbiased');

% g2=c./(mean(x))^2;
% t=[0:.000001:(length(x)-1)*0.000001];
% figure(1)
% %  raw data
% plot(t,x);
% figure(2)
% % correlation
% semilogx(lags*0.000001,c); % for laser stability
% figure(3)
% % g2 intensity autocorrelation
% semilogx(lags*.000001,g2); %for sample

% [c,lags]=xcorr(x,10000,'unbiased');


% st = time.time()
% data = pd.read_csv(filename)
% end = time.time()

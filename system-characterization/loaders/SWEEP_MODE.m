% Loads files acquired in the sweep mode:
% v.1.0
    % Inputs:
		% structure with fields:
			% "file_name": "string.csv",
			% "absorbance": float,
			% "absorbance_std": float,
			% "quantum_yield": float or "-",
			% "scattering": float,
			% "concentration": float or "-",
			% "filter_transmitance": float,
			% "solvent": {
			% 	"absorbance": float,
			% 	"name": string,
			% 	"refractive_idx": float
			% }
    % Outputs:
		% main output: quantum yield
		
	% Updates: Luminescence slope
% ========================================================================


classdef SWEEP_MODE < handle
	% SWEEP_MODE 
	%   This script loads data from the sweep mode
	properties
		absorbance
		absorbedPowAtCentre
		background % structure with fields: data, averageApd
		data
		dataDescription
		dataDetails
		details
		dataFileName
		dataFullPath
		dataPath
		expQuantumYield
		fluenceRate
		filterTransmitance
		header
		luminSignal
		powerAtTheCentre
		quantumYield
		quantumYieldNIR
		rawData
		refractiveIdx
		relQuantumYield
		scattering
		solvent % structure with fields: refractiveIdx
	end

	methods(Static)
		function [idx, v] = FindClosest(arr, value)
			% FindColsest(arr, value)
			% arr := 1D array
			% value := value to be searched
			% Returns the index and value of the closes searched number
			idx = find(abs(arr - value) == min(abs(arr - value)));
			v = arr(idx);
		end
	end 

	methods
		function DataHandler = SWEEP_MODE( paramsStruc, dataPath )
			DataHandler.dataFileName = paramsStruc.file_name;
			DataHandler.dataPath = dataPath;
			DataHandler.dataFullPath = strcat(dataPath, paramsStruc.file_name);

			DataHandler.absorbance = paramsStruc.absorbance;
			DataHandler.scattering = paramsStruc.scattering;
			DataHandler.quantumYield = paramsStruc.quantum_yield;
			DataHandler.filterTransmitance = paramsStruc.filter_transmitance;
			DataHandler.solvent.refractiveIdx = paramsStruc.solvent.refractive_idx;
			DataHandler.solvent.absorbance = paramsStruc.solvent.absorbance;

			content = importdata( DataHandler.dataFullPath );
			DataHandler.rawData = content.data;
			DataHandler.details = content.textdata;
			DataHandler.header = content.colheaders;

			% [DataHandler.rawData, ~, ~] = xlsread( DataHandler.dataFullPath );
			DataHandler.data = DataHandler.rawData(5:end-5, :);
			% Filtre current over 90 mA
			DataHandler.data = DataHandler.Filtre(110); 
			% Eliminate any row with NaN values
			DataHandler.data = DataHandler.data(any(isnan(DataHandler.data), 2)==0, :);
			
			% Get background values from initial data points
			DataHandler.background.data = DataHandler.rawData( ...
				1 : DataHandler.FindClosest(DataHandler.rawData(:, 1), 45), ...
				: ...
			);
			% Eliminate any row with NaN values
			DataHandler.background.data = DataHandler.background.data( ...
				any(isnan(DataHandler.background.data), 2) == 0, : ...
			);
			DataHandler.background.averagePower = mean(DataHandler.background.data(:, 2));
			DataHandler.background.averageApd = mean(DataHandler.background.data(:, 3));
		end 

		function data = Filtre(DataHandler, minValue)
			% Filtre() filtres the data based on a min value of the laser current
			data = DataHandler.data( DataHandler.data(:, 1) > minValue, : );
		end 
		
		function data = Data(DataHandler)
			%  Data() returns the data in a matrix
			data = DataHandler.data;
			disp(DataHandler.header);
			disp(DataHandler.data(1:20,:));
		end

		function Header(DataHandler, n)
			%  Header() returns the data in a matrix

			if ~exist('n','var')
				% optional parameter
				n = 5;
			end
			disp(DataHandler.header);
			disp(DataHandler.data(1:n,:));
		end

		function RemoveBackground(DataHandler)
			% RemoveBackground()

			DataHandler.data(:, 2) = DataHandler.data(:, 2) - DataHandler.background.averagePower;
			DataHandler.data(:, 3) = DataHandler.data(:, 3) - DataHandler.background.averageApd;
			% DataHandler.data(:, 3) = DataHandler.data(:, 3) - 0.0405;
		end

		function sample = Sample(DataHandler, nPoints)
			% Sample() returns the a sample of the original file with nPoints length
			idxs = 1;
			for i = 2 : DataHandler.Length()
				if (mod(i, round(DataHandler.Length()/(nPoints))) == 0)
					idxs = cat(2, idxs, i);
				end
			end
			sample = DataHandler.data(idxs, 1:end);
		end

		function UseSample(DataHandler, nPoints)
			% UseSample() passes a sample of the original file with nPoints length
			% to DataHandler.data
			DataHandler.data = DataHandler.Sample(nPoints);
		end

        function details = Details(DataHandler, nameDictFile)
            % Details() returns details of how the data were acquired
            % Info from dataFileName
            try
                [DataHandler.dataDescription, DataHandler.dataDetails] = nameDecoder( ...
                    DataHandler.dataFileName, ...
                    nameDictFile, ...
                    DataHandler.dataPath ...
                );
    
                details = DataHandler.dataDetails;
                disp('Data: ');
                fprintf(strcat(DataHandler.dataDescription, '\n\n'));
            catch
                warning(strcat('I couldn''t decode the file''s name.\n', ...
                    'Aparently, it isn''t written in the follwing pattern:\n', ...
                    '"ProjectName__LabsName__ResearcherInitials__YYMMDD__conditions__000.extension"\n', ...
                    'e.g.: QY-Test__T-BioP2__JM__180915__sTest-apd1-g_8__001.png' ...
                ));
            end
        end

		function powerAtTheCentre = PowAtCentre(DataHandler)
			% 
			DataHandler.powerAtTheCentre = DataHandler.LaserPower() ...
				.* sqrt(exp(DataHandler.absorbance + DataHandler.solvent.absorbance));
			powerAtTheCentre = DataHandler.powerAtTheCentre;
		end

		function absorbance = Absorbance(DataHandler, emptyCuvetteHolderObj, returnType) 
			% Absorbance(laserBehaviour) returns absorbance of the sample
			% absorbance = log(power before cuvette / power after cuvette)
			% NOT IN USE

			laserPower = polyval( ...
				polyfit( ...
					emptyCuvetteHolderObj.LaserCurrent(), ...
					emptyCuvetteHolderObj.LaserPower(), ...
					1 ...
				), ...
				DataHandler.LaserCurrent() ...
			);
			if returnType == 'mean' 
				absorbance = log(laserPower ./ DataHandler.LaserPower());
				absorbance = mean(absorbance (absorbance > 0));
			elseif returnType == 'vector'
				absorbance = log(laserPower ./ DataHandler.LaserPower());
			end
		end

		function LuminSignal(DataHandler) 
            % Luminescent signal of a suspension
            DataHandler.luminSignal = DataHandler.APD1Voltage();
		end

		function GetLuminSlope(DataHandler)
			% Calculates the slope of the luminescence curve in respect to absorbed power
			x = log10(DataHandler.fluenceRate);
			y = log10(DataHandler.luminSignal);
			% % y = rand (lengthX,1);
			% % Plot it and show how the line has sharp bends.
			samplingRateIncrease = 1;
			newXSamplePoints = linspace(x(1), x(end), length(x) * samplingRateIncrease);
			smoothedY = spline(x, log10(DataHandler.luminSignal), newXSamplePoints);
			slopes = [0, diff(smoothedY)./diff(newXSamplePoints)];
			smoothSlopes = movmean(slopes, 35);
			% % slopes = [0, diff(log(DataHandler.luminSignal)) ./ diff(log(DataHandler.fluenceRate))];
			% % Plot smoothedY and show how the line is
			
			% dydx = (y(2:end) - y(1:end-1)) ./ (x(2:end) - x(1:end-1));
			% x = x(2:end);
			figure(2)
				subplot(2, 1, 1)
					plot(log10(DataHandler.fluenceRate), log10(DataHandler.luminSignal), '-sr', 'LineWidth', 2);
					hold on; % Don't destroy the first curve we plotted.
					plot(newXSamplePoints, smoothedY, '-ob');
					ylabel('log(I)')
					xlabel('log(Power density)')
					title('Slope of log(I) x log(power density)', 'FontSize', 20);
					legend('Original Points', 'Spline Points');
			
				subplot(2, 1, 2)
					plot(newXSamplePoints(10:end), slopes(10:end), 'ko-', 'MarkerSize', 3);
					hold on; % Don't destroy the first curve we plotted.
					plot(newXSamplePoints(10:end), smoothSlopes(10:end), 'r.-', 'MarkerSize', 3);
					xlabel('log(Power density)')
					ylabel('d/d[log(p)] x log(I)')
					% plot(x, dydx, 'gd-', 'MarkerSize', 3);
					% Draw x axis
					% line('Color', 'k', 'LineWidth', 2);
					grid on;
					legend('Slope', 'Moving average 35');
		end

		function LuminSlope(DataHandler)
			% 
		end

		function absorbedPowAtCentre = AbsorbedPowAtCentre(DataHandler) 
			% AbsorbedPowAtCentre(emptySampleHolderObj) multiplies absorbance by power intensity 
			% at the centre of the cuvette. 

			DataHandler.absorbedPowAtCentre = (DataHandler.absorbance - DataHandler.scattering) ...
				.* DataHandler.PowAtCentre(); 
			absorbedPowAtCentre = DataHandler.absorbedPowAtCentre;
		end

		function expQuantumYield = ExpQuantumYield(DataHandler) 
			% 

			DataHandler.expQuantumYield = DataHandler.luminSignal ./ DataHandler.absorbedPowAtCentre;
			expQuantumYield = DataHandler.expQuantumYield;
		end

		function relQuantumYield = RelQuantumYield(DataHandler, REFERENCE)
			% RelQuantumYield - Calculates the relative quantum yield
			%
			% Syntax: relQuantumYield = RelQuantumYield(REFERENCE)
			% .8387 -> filter transmitance of UCNPs
			% .35 -> filter transmitance of Dye-781
			% .124 -> QY of dye
			qyRef = polyfit(REFERENCE.absorbedPowAtCentre, REFERENCE.luminSignal, 1);
			disp(qyRef)
			DataHandler.relQuantumYield = ( ...
				((DataHandler.ExpQuantumYield() / DataHandler.filterTransmitance) * REFERENCE.quantumYield) ./ ... 
				(qyRef(1) / REFERENCE.filterTransmitance) ...
			) * ((DataHandler.solvent.refractiveIdx / REFERENCE.solvent.refractiveIdx) ^ 2);

			relQuantumYield = DataHandler.relQuantumYield;
		end

		function fittingParams = QYFitting(DataHandler, model, guessArray)
			%QYFitting - Fits QY experimental data with model
			%
			% Syntax: fittingParams = QYFitting(model, guessArray)
			%
			% Long description

			options = optimoptions( ...
				'lsqcurvefit', 'OptimalityTolerance', 1e-20, ...
				'FunctionTolerance', 1e-20, 'StepTolerance', 1e-10, ...
				'MaxFunctionEvaluations', 500 ...
			);

			fittingParams = lsqcurvefit( ...
				model, ...
				guessArray, ...
				DataHandler.fluenceRate, ... 
				DataHandler.relQuantumYield, ...
				[], ...
				[], ...
				options ...
			);
			disp('QY Saturation: ');
			fprintf(strcat(num2str(fittingParams(1)), '\n\n'));
			disp('QY Balancing Point(s): ');
			fprintf(strcat(num2str(fittingParams(2:end)), '\n\n'));
		end

		function fluenceRate = FluenceRate(DataHandler, beamArea) 
			% 

			DataHandler.fluenceRate = ( ...
				((1 - 0.5 * (DataHandler.absorbance - DataHandler.scattering)) ... 
				.* DataHandler.powerAtTheCentre) / beamArea ...
			);

			fluenceRate = DataHandler.fluenceRate;
		end

		function col1 = LaserCurrent(DataHandler)
			col1 = DataHandler.data(:, 1);
		end

		function col2 = LaserPower(DataHandler)
			col2 = DataHandler.data(:, 2) ./ 1000; % convert to mW
		end

		function col3 = APD1Voltage(DataHandler)
			col3 = DataHandler.data(:, 3);
		end

		function len = Length(DataHandler)
			% Length returns the length of the 
			% dataset based on the APD1 Voltage column
			len = size(DataHandler.APD1Voltage(), 1);			
		end

		function Plot(DataHandler, which, legendLabel)
			% Plot(which, legendLabel)
			% legendLabel is optional 
			% which = [rawdata, expQY, relQY, luminescence]

			if ~exist('legendLabel','var')
				% optional parameter
				legendLabel = split( ...
					DataHandler.dataFileName, ...
					'__' ...
				);
				legendLabel = legendLabel{1};
			end

			% colours = {'b', 'm', 'r', 'c', 'k', 'g', 'y'};
			shapes = {'^', 'o', 'x', 's', '*', 'd', '+', 'v', 'p', '>', 'h', '<'};

			% colourShape = strcat(colours{randi(numel(colours))}, shapes{randi(numel(shapes))});
			shape = shapes{randi(numel(shapes))};

			if strcmp(which, 'rawdata')
				f1 = figure(100);
				screenSize = get(0, 'ScreenSize');
				set(f1, 'position', [900, 50, screenSize(3) / 3.4, screenSize(4) / 1.15]);
					subplot(3, 1, 1);
					plot( ...
						DataHandler.rawData(1:end-1, 1), ...
						DataHandler.rawData(1:end-1, 3), ...
						'-o', ...
						'DisplayName', legendLabel ...
					);
					% DataHandler.LaserCurrent(), ...
					% DataHandler.APD1Voltage(), ...
					hold on
					% xlim([50,max(DataHandler.rawData(1:end-1, 1))])
					% title(sampleName)
					xlabel('I (mA)');
					ylabel('V (mV)');
					legend('location', 'best');
				
					subplot(3, 1, 2)
					plot( ...
						DataHandler.rawData(1:end-1, 1), ...
						DataHandler.rawData(1:end-1, 2), ...
						'-s', ...
						'DisplayName', legendLabel ...
					);
					
						% DataHandler.LaserCurrent, ...
						% DataHandler.LaserPower, ...
					hold on
					% xlim([50, max(DataHandler.LaserCurrent)])
					xlabel('I (mA)');
					ylabel('P (uW)');
				
					subplot(3, 1, 3);
					plot( ...
						DataHandler.rawData(1:end-1, 2), ...
						DataHandler.rawData(1:end-1, 3), ...
						'-d', ...
						'DisplayName', legendLabel ...
					);
						% DataHandler.LaserPower, ...
						% DataHandler.APD1Voltage, ...
					
					hold on
					% xlim([0, max(DataHandler.LaserPower)])
					xlabel('P (uW)');
					ylabel('V (mV)');

			elseif strcmp(which, 'expQY')
				f2 = figure(102);
					plot( ...
						DataHandler.fluenceRate, ...
						DataHandler.expQuantumYield, ...
						shape, ...
						'DisplayName', legendLabel ...
					);
					
					hold on
					title('Exprimental QY')
					xlabel('Average power density (mW/cm^2)')
        			ylabel('Quantum Yield (-)')
					legend('location', 'best');

			elseif strcmp(which, 'relQY')
				f3 = figure(103);
					plot( ...
						DataHandler.fluenceRate, ...
						DataHandler.relQuantumYield, ...
						shape, ...
						'MarkerSize', 5, ...
						'DisplayName', legendLabel ...
					);
					
					hold on
					title('Relative QY')
					xlabel('Average power density (mW/cm^2)')
					ylabel('Quantum Yield (-)')
					legend('location', 'best');

			elseif strcmp(which, 'luminescence')
				figure(104)
					plot( ...
						DataHandler.absorbedPowAtCentre, ...
						DataHandler.luminSignal, ...
						shape, ...
						'MarkerSize', 5, ...
						'DisplayName', legendLabel ...
					);
					
					hold on
					title('Luminescence signal')
					xlabel('Absorbed power(mW/cm)')
					ylabel('Luminescence signal [counts]')
					legend('show', 'location', 'best');
			end
		end
	end
end
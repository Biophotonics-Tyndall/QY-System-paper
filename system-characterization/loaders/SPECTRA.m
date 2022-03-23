% Quantum Yield - Emission over time loader class:
    % Inputs:
        % data file full path
    % Outputs:
        %
% ========================================================================


classdef SPECTRA < handle
	% SPECTRA 
	%   This script loads data from the sweep mode
	properties
		absorptionSpectrum
		BACKGROUND
		cuvetteLength
		data
		dataDescription
		dataDetails
		dataFileName
		dataFullPath
		dataPath
		details
		laserCurrent
		rawData
		REFERENCE
	end

	methods(Static)
		function idx = FindClosest(arr, value)
			idx = find(abs(arr - value) == min(abs(arr - value)));
		end
	end 

	methods
		function DataHandler = SPECTRA( dataFileName, dataPath )
			DataHandler.dataFileName = dataFileName;
			DataHandler.dataPath = dataPath;
			DataHandler.dataFullPath = strcat(dataPath, dataFileName);

			content = importdata(DataHandler.dataFullPath);
			DataHandler.rawData = content.data;
			DataHandler.data = DataHandler.rawData;
			DataHandler.details = content.textdata;
			DataHandler.cuvetteLength = 1; % cm
		end
		
		function data = Data(DataHandler)
			%  Data() returns the data in a matrix
			data = DataHandler.data;
		end

		function data = SetReference(DataHandler, REFERENCE)
			%SetReference - Description
			%
			% Syntax: SetReference(DataHandler, REFERENCE)
			%
			% Long description
			DataHandler.REFERENCE = REFERENCE;
			data = REFERENCE.Intensity() - DataHandler.data(:, 2);
			DataHandler.data(:, 3) = data;
		end

		function data = SetBackground(DataHandler, BACKGROUND)
			%SetBackground - Description
			%
			% Syntax: SetBackground(DataHandler, BACKGROUND)
			
			DataHandler.BACKGROUND = BACKGROUND;
			data = DataHandler.data(:, 2) - BACKGROUND.Intensity(); 
			DataHandler.data(:, 2) = data;
		end
		
		function mu_a = GetAbsorbance(DataHandler)
			%GetAbsorbance - Description
			%
			% Syntax: data = GetAbsorbance(DataHandler)
			%
			% reference := mu_a
			mu_a = log( ...
				DataHandler.REFERENCE.Intensity() ./ DataHandler.Intensity() ...
			) ./ DataHandler.cuvetteLength;
			DataHandler.data(:, 4) = mu_a;
		end

		function [abArr, model] = FitAbsorbance(DataHandler, params)
			% FitAbsorbance - Fit absorbance with MODELS.Absorbance method
			% taking params as initial guesses and parameters for the chosen model
			% returns array abArr := [mu_0, b] = FitAbsorbance(DataHandler, params)
			
			from = DataHandler.FindClosest(DataHandler.data(:, 1), params.from);
			to = DataHandler.FindClosest(DataHandler.data(:, 1), params.to);

			model = @(abArr, lambda) MODELS.Absorbance(abArr, lambda, params.lambda0);

            options = optimoptions( ...
				'lsqcurvefit', 'OptimalityTolerance', 1e-20, ...
				'FunctionTolerance', 1e-20, 'StepTolerance', 1e-10, ...
				'MaxFunctionEvaluations', 500 ...
            );
            
			[abArr, residuals] = lsqcurvefit( ...
				model, ...
				[params.initial_guess.mu_0, params.initial_guess.b], ...
				DataHandler.data(from :to, 1), ... 
				DataHandler.data(from :to, 4), ...
				[], ...
                [], ...
                options ...
			);
			disp(abArr);
			disp(residuals);

			wavelength = linspace( ...
				DataHandler.data(10,1), ...
				DataHandler.data(end-10,1), ...
				300 ...
			).';

			fittedAbsorbance = model(abArr, wavelength);
			from = DataHandler.FindClosest(wavelength, params.from);
			to = DataHandler.FindClosest(wavelength, params.to);

			figure(202)
			plot( ...
				wavelength, ...
				fittedAbsorbance, ...
				'r-', ...
				'DisplayName', ...
				'Fitting' ...
			);
			hold on
			plot( ...
				wavelength(from :to), ...
				fittedAbsorbance(from :to), ...
                'c-', ...
                'DisplayName', ...
				'Region used' ...
			);
			
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

		function Normalize(DataHandler) 
			% 
			DataHandler.data(:, 2) = DataHandler.data(:, 2) ./ max(DataHandler.data(end - 250 : end, 2));
		end

		function NormalizeMod(DataHandler, arrRange) 
			% 
			if ~exist('arrRange','var')
				% optional parameter
				DataHandler.data(:, 2) = DataHandler.data(:, 2) ./ max(DataHandler.data(:, 2));
			else
				from = DataHandler.FindClosest(DataHandler.Wavelength(), arrRange(1));
				to = DataHandler.FindClosest(DataHandler.Wavelength(), arrRange(2));
				DataHandler.data(:, 2) = DataHandler.data(:, 2) ./ max(DataHandler.data(from: to, 2));
			end
		end

		function RemoveBaseLine(DataHandler) 
			%
			% DataHandler.data(:, 2) = DataHandler.data(:, 2) - min(DataHandler.data(20:end-50, 2));
			baseLineData = cat( ...
				1, ...
				DataHandler.data(20: 90, :), ...
				DataHandler.data(end-140: end-50, :) ...
			);
			baseLine = polyval( ...
				polyfit( ...
					baseLineData(:, 1), ...
					baseLineData(:, 2), ...
					1 ...
				), ...
				DataHandler.Wavelength() ...
			);

			DataHandler.data(:, 2) = DataHandler.data(:, 2) - baseLine;
			DataHandler.data(:, 2) = DataHandler.data(:, 2) - min(DataHandler.data(20:end-50, 2)) + 0.001;

		end

		
		function RemoveBaseLineMod(DataHandler, arrayRange) 
			%
			% DataHandler.data(:, 2) = DataHandler.data(:, 2) - min(DataHandler.data(20:end-50, 2));

			from1 = DataHandler.FindClosest(DataHandler.Wavelength(), arrayRange(1));
			to1 = DataHandler.FindClosest(DataHandler.Wavelength(), arrayRange(2));
			from2 = DataHandler.FindClosest(DataHandler.Wavelength(), arrayRange(3));
			to2 = DataHandler.FindClosest(DataHandler.Wavelength(), arrayRange(4));
			
			baseLineData = cat( ...
				1, ...
				DataHandler.data( from1: to1, :), ...
				DataHandler.data(from2: to2, :) ...
			);
			baseLine = polyval( ...
				polyfit( ...
					baseLineData(:, 1), ...
					baseLineData(:, 2), ...
					1 ...
				), ...
				DataHandler.Wavelength() ...
			);

			DataHandler.data(:, 2) = DataHandler.data(:, 2) - baseLine;
			% DataHandler.data(:, 2) = DataHandler.data(:, 2) - min(DataHandler.data(20:end-50, 2)) + 0.001;

		end

		function col1 = Wavelength(DataHandler)
			col1 = DataHandler.data(:, 1);
		end

		function col2 = Intensity(DataHandler)
			col2 = DataHandler.data(:, 2);
		end

		function col3 = AbsorbedPow(DataHandler)
			col3 = DataHandler.data(:, 3);
		end

		function col4 = Absorbance(DataHandler)
			col4 = DataHandler.data(:, 4);
		end

		function col5 = FittedAbsorbance(DataHandler)
			col5 = DataHandler.data(:, 5);
		end

		function len = Length(DataHandler)
			% Length returns the length of the 
			% dataset based on the APD1 Voltage column
			len = size(DataHandler.Wavelength(), 1);			
		end



		function Plot(DataHandler, which, legendLabel)

			if ~exist('legendLabel','var')
				% optional parameter
				legendLabel = split( ...
					DataHandler.dataFileName, ...
					'__' ...
				);
				legendLabel = legendLabel{1};
			end

			if strcmp(which, 'intensity')
				figNum = 200;
				y = DataHandler.Intensity();
				yLabel = 'Intensity (counts)';

			elseif strcmp(which, 'absorbedPow')
				figNum = 201;
				y = DataHandler.AbsorbedPow();
				yLabel = 'Absorbed Power (a.u.)';
				
			elseif strcmp(which, 'absorbance')
				figNum = 202;
				y = DataHandler.Absorbance();
				yLabel = 'Absorbance (a.u.)';
			end

			figure(figNum)
				plot( ...
					DataHandler.Wavelength(), ...
					y, ...
					'.', ... 
					'DisplayName', ...
					legendLabel ...
				);
				hold on
				xlim([350,1100])
				% title(sampleName)
				xlabel('Wavelength (nm)');
				ylabel(yLabel);
				% legend;
		end
	end
end
% Quantum Yield - Emission over time loader class:
    % Inputs: 
        % file name
    % Outputs:
        % matrix of pixel intensities
% ========================================================================

classdef MODELS %< handle
    % IntensityVsTime Summary of this class goes here
    % This script loads and remove background from data
    % properties
    %     data
    % end
    
    methods(Static)   
        function fittingParams = QYFitting(UCNP, model, guessArray)
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
				UCNP.fluenceRate, ... 
				UCNP.relQuantumYield, ...
				[], ...
				[], ...
				options ...
			);
			disp('QY Saturation: ');
			fprintf(strcat(num2str(fittingParams(1)), '\n\n'));
			disp('QY Balancing Point(s): ');
            fprintf(strcat(num2str(fittingParams(2:end)), '\n\n'));
        end
        
        function eta = QuantumYieldNIR(abArr, powerDens, normalizedBeamProfile, beamSpotArea, pxArea)
            % Quantum yield for a beam profile distribution given by normalizedBeamProfile
            % normalizedBeamProfile has total volume = 1. 
            % It means that the its intensity has arbitrary value 
            
            [n, ~] = size(powerDens);
            eta = zeros(n, 1); % Initializing the uniform power density quantum yield vector
            powerDens = powerDens .* beamSpotArea ./ pxArea;
            for k = 1 : n      
                eta(k) = powerDens(k) .* abArr(1) * sum( ...
                    sum( ...
                        ( ...
                            normalizedBeamProfile .* normalizedBeamProfile ./ ...
                            (abArr(2) + powerDens(k) .* normalizedBeamProfile) ...
                        ) ...
                    ) ...
                );
            end
        end

        function mu_a = Absorbance(abArr, lambda, lambda_0)
            % abArr = [mu_0, b]
            mu_a = abArr(1) .* (lambda ./ lambda_0) .^ (-abArr(2));
        end

    end
end
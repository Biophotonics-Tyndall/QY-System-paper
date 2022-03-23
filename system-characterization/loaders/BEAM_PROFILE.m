% Quantum Yield - Emission over time loader class:
    % Inputs: 
        % file name, name of the dictionary file to decode the data file name, ...
        % path to data, camera pixel size (usually in cm)
    % Figures start at 200
% ========================================================================

classdef BEAM_PROFILE < handle
    % BEAM_PROFILE is class to handle laser beams profiles
        % Each function bellow (apart from BEAM_PROFILE()) represents a method
        % to analyse the beam profile of a laser spot. Each method must be 
        % documented individually.

    properties % are initialized with the constructor method
        % provided by user:
        dataFileName 
        dataPath % path where data is stored. (recommended to keep them on ../data/)
        camPxSize % camera's pixel size provided by user
        beamWidth
        widthMethod

        % evaluate from info provided
        dataFullPath % is: dataPath/dataFileName
        data % is modified everytime a new intended method is called
        rawData % original data without any modification
        dataDescription % info from decodification of the dataFileName
        dataDetails % keeps the information from the filename after being decoded

    end
    
    methods
        function DataHandler = BEAM_PROFILE(dataFileName, dataPath, widthMethod, camPxSize)
            % BEAM_PROFILE is the constructor method, i.e. it is run automatically 
            % when the class is called for the first time in an external script.
            % In other words, it initializes the the data object. 
            % For more info on classes, check this link out: https://bit.ly/2Ez3GTN
            
            DataHandler.widthMethod = widthMethod;
            DataHandler.camPxSize = camPxSize;

            DataHandler.dataFileName = dataFileName;
            DataHandler.dataPath = dataPath;
            DataHandler.dataFullPath = strcat (dataPath, dataFileName);
            DataHandler.rawData = imread( DataHandler.dataFullPath );
            try
                DataHandler.data = double(rgb2gray(DataHandler.rawData)); 
            catch
                DataHandler.data = double(DataHandler.rawData);
            end
           
        end
        
        function data = Data(DataHandler)
            % Data() returns matrix of pixel intensities

            data = DataHandler.data;
        end
        
        function Reset(DataHandler) 
            % Reset() method replace DataHandler.data by DataHandler.rawData

            try
                DataHandler.data = double(rgb2gray(DataHandler.rawData)); 
            catch
                DataHandler.data = double(DataHandler.rawData);
            end
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

        function RemoveBackground(DataHandler) 
            % RemoveBackground() finds a region away from the beam spot and  
            % uses it as a sample of background noise. Its max value is 
            % used as a threshold to set the pixes' intensities to 0.

            matrixSize = size(DataHandler.data);
            beamWidthXY = DataHandler.Width(1);
            if matrixSize(1) < 4 * beamWidthXY || matrixSize(2) < 4 * beamWidthXY
                disp('Warning: Background not removed! Did you trim the data matrix?')
            else
                % finds limits of the noise sample
                spotCentre = round(DataHandler.MaxPeak()); % 
                [rowLength, colLength] = size(DataHandler.data); 
                rowMax = spotCentre(2) + 3 * beamWidthXY;
                rowMin = spotCentre(2) + 2 * beamWidthXY;
                colMax = spotCentre(3) + 3 * beamWidthXY;
                colMin = spotCentre(3) + 2 * beamWidthXY;
                if rowMax >= rowLength 
                    rowMin = spotCentre(2) - 3 * beamWidthXY;
                    rowMax = spotCentre(2) - 2 * beamWidthXY;
                end
                if colMax >= colLength 
                    colMin = spotCentre(3) - 3 * beamWidthXY;
                    colMax = spotCentre(3) - 2 * beamWidthXY;
                end
    
                backgroundLevel = max(max(DataHandler.data(rowMin : rowMax, colMin : colMax)));
                DataHandler.data( DataHandler.data <= backgroundLevel) = 0;
            end
        end

        function Trim(DataHandler, lengthScale)
            % Trim(lengthScale) trims the image around the laser spot
            % Set lengthScale from 0.5 up to 2. This scales the 
            % new matrix proportionally to the beam width

            spotCentre = DataHandler.MaxPeak();
            increment = lengthScale * DataHandler.Width(1); 
            rowMax = round(spotCentre(2) + increment);
            rowMin = round(spotCentre(2) - increment);
            colMax = round(spotCentre(3) + increment);
            colMin = round(spotCentre(3) - increment);
            DataHandler.data = DataHandler.data(rowMin : rowMax, colMin : colMax);
        end

        function Normalize(DataHandler) 
            % Normalize() method normalizes the matrix in order the volume 
            % under the curve is 1

            DataHandler.data = DataHandler.data ./ sum(sum(DataHandler.data));
        end

        function maxPeakArr = MaxPeak(DataHandler)
            % MaxPeak() returns a coord array with the most intense px
            % [maxValue, row, col]
            
            [maxNumCol, maxIndexCol] = max(DataHandler.data);
            [maxNum, col] = max(maxNumCol);
            row = maxIndexCol(col);
            maxPeakArr = [double(maxNum), row, col];
        end
        
        function [beamWidthY, beamWidthX] = Width(DataHandler, pxSize) 
            % Width() finds the most bright pixel and returns an array for x and y widths
            % If pxSize is passed as 1, it returns the size in number of pixels.
            % Two methods can be selected: '1/e2' or 'fwhm'
            
            if DataHandler.widthMethod == '1/e2' 
                coeff = 0.135; 
            elseif DataHandler.widthMethod == 'fwhm'
                coeff = 0.5; 
            end
            maxPeakArr = DataHandler.MaxPeak();
            row = maxPeakArr(2);
            col = maxPeakArr(3);
            maxNum = maxPeakArr(1);
            coeff = coeff * maxNum;
            % vertical
            [~, widthY] = size(find(DataHandler.data(row, :) > coeff));
            beamWidthY = widthY * pxSize;
            % horizontal
            [widthX, ~] = size(find(DataHandler.data(:, col) > coeff));
            beamWidthX = widthX * pxSize;

            DataHandler.beamWidth = (beamWidthX + beamWidthY) ./ 2;
        end

        function spotArea = Area(DataHandler) 
            % Area() returns the area of the beam spot
            % widthMethods: 1/e2 or fwhm

            spotArea = pi * (mean(DataHandler.Width(DataHandler.camPxSize)) / 2) ^ 2;
        end
        
        function gaussMatrix = GaussianMatrix(DataHandler)
            % GaussianMatrix() returns a normalized gaussian matrix with 
            % experimental beam width 
            % width methods: '1/e2' or 

            beamWidthXY = DataHandler.Width(DataHandler.camPxSize);
            [m, l] = size(DataHandler.data);
            squareSideLength = DataHandler.camPxSize * l; % size of the matrix (It can be reduced to the size of the beam spot)
            mVec = linspace(-0.5 * squareSideLength, 0.5 * squareSideLength, l);  
            gaussMatrix = zeros([m l]);
            % Defining a normalized beam profile based on the beam profiler image 
            % This for loop can be removed! 
            for m_i = 1 : m
                for l_i = 1 : l
                    % Defining a basic Gaussian profile with the experimental beam width
                    gaussMatrix(m_i, l_i) = exp( ...
                    -(1 / sqrt(2)) * ...
                    (mVec(m_i) ^ 2 + mVec(l_i) ^ 2) / ...
                    (0.25 * beamWidthXY ^ 2) ...
                    );             
                end
            end
            gaussMatrix = gaussMatrix ./ sum(sum(gaussMatrix));
        end
        
        function Imshow(DataHandler)
            % ImShow() displays the original picture, i.e. from DataHandler.rawData

            figure(200)
            imshow(DataHandler.rawData)
        end
        
        function Surf(DataHandler, m)
            % Surf plots a 3D visualization of the beam profile
            
            
            if m == "image" 
                centre = double(DataHandler.MaxPeak());
                beamWidthXY = DataHandler.Width(1); % in pxs
                xLimArr = [centre(3) - 1.25 * beamWidthXY(1), centre(3) + 1.25 * beamWidthXY(1)];
                yLimArr = [centre(2) - 1.25 * beamWidthXY(1), centre(2) + 1.25 * beamWidthXY(1)];
                
                figure (201);
                s = surf(DataHandler.data);
                xlim(xLimArr);
                ylim(yLimArr);
            elseif m == "gaussian"
                figure (202);
                s = surf(DataHandler.GaussianMatrix());
            end
            set(s, 'edgecolor','none');
            set(get(gca, 'XLabel'), 'String', '[px]');
            set(get(gca, 'YLabel'), 'String', '[px]');
            set(get(gca, 'ZLabel'), 'String', 'Brightness [counts]');
        end
    end
end







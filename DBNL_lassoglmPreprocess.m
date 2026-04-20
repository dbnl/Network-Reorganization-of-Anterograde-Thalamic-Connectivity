function [processedMatrix, variableNames, responseVar, nanRows] = DBNL_lassoglmPreprocess(dataTable, responseColumn)
    % DBNL_lassoglmPreprocess - Dynamically preprocess a dataset for lassoglm
    %
    % Inputs:
    %   dataTable - Input table with predictors and response variable
    %   responseColumn - Name of the response variable column (string or char)
    %
    % Outputs:
    %   processedMatrix - Numeric matrix of preprocessed predictors
    %   variableNames - Cell array of predictor names
    %   responseVar - Vector of response variable values (numeric, for lassoglm)
    %   nanRows - Logical vector indicating rows removed due to NaNs

    % Initialize storage for processed data
    processedMatrix = [];
    variableNames = {};

    % Separate response variable
    responseVar = dataTable.(responseColumn);
    if islogical(responseVar) || iscategorical(responseVar)
        responseVar = double(responseVar); % Convert to numeric if binary
    end

    % Iterate through each predictor
    predictorColumns = setdiff(dataTable.Properties.VariableNames, responseColumn);
    for col = predictorColumns
        predictorData = dataTable.(col{1}); % Extract column data
        
        % Handle different data types dynamically
        if isnumeric(predictorData) || islogical(predictorData)
            % Handle numeric or logical data
            if all(ismember(unique(predictorData), [0, 1])) % Binary
                processedMatrix = [processedMatrix, double(predictorData)]; % Convert logical to double
                variableNames = [variableNames, col{1}];
            else
                % Normalize continuous data
                predictorMean = mean(predictorData, 'omitnan');
                predictorStd = std(predictorData, 'omitnan');
                if predictorStd == 0
                    normalizedData = zeros(size(predictorData)); % Handle constant columns
                else
                    normalizedData = (predictorData - predictorMean) / predictorStd;
                end
                processedMatrix = [processedMatrix, normalizedData];
                variableNames = [variableNames, col{1}];
            end
        elseif iscellstr(predictorData) || isstring(predictorData) || iscategorical(predictorData)
            % Handle categorical data
            catData = categorical(predictorData); % Ensure it's categorical

            % Group rare categories into 'other'
            counts = countcats(catData);
            total = sum(counts);
            rareCategories = categories(catData);
            rareCategories = rareCategories(counts / total < 0.01);
            if ~isempty(rareCategories)
                catData = mergecats(catData, rareCategories, 'other');
            end

            % Generate dummy variables
            dummyData = dummyvar(catData);
            
            % Append dummy variables and names
            processedMatrix = [processedMatrix, dummyData];
            dummyVarNames = strcat(col{1}, '_', string(categories(catData)));
            variableNames = [variableNames, dummyVarNames'];
        elseif iscell(predictorData)
            % Ensure cell array contains only character vectors
            if all(cellfun(@ischar, predictorData)) || all(cellfun(@isstring, predictorData))
                predictorData = categorical(string(predictorData)); % Convert to string then to categorical
            else
                warning('Skipping column "%s": contains non-character data.', col{1});
                continue;
            end
            
            % Group rare categories into 'other'
            counts = countcats(predictorData);
            total = sum(counts);
            rareCategories = categories(predictorData);
            rareCategories = rareCategories(counts / total < 0.01);
            if ~isempty(rareCategories)
                predictorData = mergecats(predictorData, rareCategories, 'other');
            end
            
            % Convert categorical to dummy variables
            dummyData = dummyvar(predictorData);
        
            % Append dummy variables and names
            processedMatrix = [processedMatrix, dummyData];
            dummyVarNames = strcat(col{1}, '_', string(categories(predictorData)));
            variableNames = [variableNames, dummyVarNames'];
        else
            % Unsupported type: Warn and skip the column
            warning('Skipping unsupported data type for column: %s', col{1});
        end
    end
    % Remove rows with NaN values
    nanRows = any(isnan([processedMatrix, responseVar]), 2); % Identify rows with NaNs
    processedMatrix(nanRows, :) = []; % Remove NaN rows from predictors
    responseVar(nanRows) = []; % Remove corresponding rows from response variable
end


%% Train KNN Model ( supervised learning)
% Load the spikes for brightest and dimmest 
% Load the data
Data = load('/Users/Ines/Dropbox/Retinotectal_analysis_for Nilou/matlab_analysis/R22-27_analysis_2023/C1_YASS/workspace.mat');
avg_D_cell=load('spike_average_epoch1.mat')
avg_B_cell=load('spike_average_epoch9.mat')
Algorithm_identified_ReaChR= load('ReachR_cells_computationally.mat');
Algorithm_identified_NonReaChR= load('Non_ReachR_cells_computationally.mat');

% load spike times for all the cells 
spike_times = Data.data{1,1}.spikes;

% number of cells
num_cells = length(spike_times);

% load the trigger times 
trigger_times = Data.data{1,1}.triggers;

% Define the flash duration sets 
flash_durations = Data.flash_dur;

% Start and End trigger times for each epoch (s) 
epoch_times = Data.trigtimes;

%cell IDs
cell_ids = Data.data{1,1}.cell_ids;
cells=1:length(cell_ids);

% Clustering data
avg_D_cell=avg_D_cell.avg_D_cell;
avg_B_cell=avg_B_cell.avg_B_cell;

Algorithm_identified_ReaChR=Algorithm_identified_ReaChR.ReachR_cells;
Algorithm_identified_NonReaChR=Algorithm_identified_NonReaChR.Non_ReachR_cells;


clustering_data=[avg_D_cell', avg_B_cell'];
k=3;

% Find the cell numbers for Algorithm_identified_ReaChR
[~, cell_numbers_ReaChR] = ismember(Algorithm_identified_ReaChR, cell_ids);

% Find the cell numbers for Algorithm_identified_NonReaChR
[~, cell_numbers_NonReaChR] = ismember(Algorithm_identified_NonReaChR, cell_ids);


% Create label vectors
% Initialize a label vector with zeros for all cell numbers
labels = zeros(length(cells), 1);

% Assign label 1 for ReaChR cells
labels(cell_numbers_ReaChR) = 1;

% Assign label 0 for Non-ReaChR cells
labels(cell_numbers_NonReaChR) = 2;



rng('default');  % For reproducibility
cv = cvpartition(labels, 'HoldOut', 0.3);
idx = cv.test;

% Separate the data into training and testing sets
dataTrain = clustering_data(~idx, :);
dataTest = clustering_data(idx, :);
labelsTrain = labels(~idx);
labelsTest = labels(idx);


% Create the KNN classifier
knnModel = fitcknn(dataTrain, labelsTrain, 'NumNeighbors', k);

% Make predictions on the test set
predictedLabels = predict(knnModel, dataTest);

% Calculate the classification accuracy
accuracy = sum(predictedLabels == labelsTest) / numel(labelsTest);

% Display the accuracy
fprintf('KNN classification accuracy: %.2f%%\n', accuracy * 100);

% Get the range of the data
x_min = min(clustering_data(:, 1));
x_max = max(clustering_data(:, 1));
y_min = min(clustering_data(:, 2));
y_max = max(clustering_data(:, 2));

% Create a grid of points
[x, y] = meshgrid(linspace(x_min, x_max, 100), linspace(y_min, y_max, 100));
grid_points = [x(:), y(:)];

predicted_grid_labels = predict(knnModel, grid_points);

% Plot the decision boundaries
figure();
gscatter(grid_points(:, 1), grid_points(:, 2), predicted_grid_labels, 'rb', '.', 10);
hold on;

% Plot the data points
gscatter(dataTrain(:, 1), dataTrain(:, 2), labelsTrain, 'rb', 'o', 10);

% Customize the plot
xlabel('Feature 1');
ylabel('Feature 2');
title('KNN Decision Boundaries and Data Points');
legend('ReaChR', 'Non-ReaChR', 'Location', 'best');
hold off;



save('knnModel.mat', 'knnModel');





%%  find the optimal value of K
% Choose the number of nearest neighbors (e.g., 3)
k_values = 1:20;

% Set the number of folds for cross-validation
num_folds = 5;


% Initialize a matrix to store the accuracy for each fold and K-value
accuracy_matrix = zeros(length(k_values), num_folds);


% Perform K-fold cross-validation
cv = cvpartition(labels, 'KFold', num_folds);



for k_idx = 1:length(k_values)
    k = k_values(k_idx);

    for fold_idx = 1:num_folds
        % Get the training and validation sets for the current fold
        train_idx = cv.training(fold_idx);
        validation_idx = cv.test(fold_idx);
        dataTrain = clustering_data(train_idx, :);
        dataValidation = clustering_data(validation_idx, :);
        labelsTrain = labels(train_idx);
        labelsValidation = labels(validation_idx);

        % Train the KNN classifier
        knnModel = fitcknn(dataTrain, labelsTrain, 'NumNeighbors', k);

        % Make predictions on the validation set
        predictedLabels = predict(knnModel, dataValidation);

        % Calculate the classification accuracy
        accuracy = sum(predictedLabels == labelsValidation) / numel(labelsValidation);
        
        % Store the accuracy in the accuracy_matrix
        accuracy_matrix(k_idx, fold_idx) = accuracy;
    end
end

% Calculate the average accuracy for each K-value
mean_accuracy = mean(accuracy_matrix, 2);

% Find the optimal K-value
[~, optimal_k_idx] = max(mean_accuracy);
optimal_k = k_values(optimal_k_idx);

fprintf('The optimal K-value is: %d\n', optimal_k);


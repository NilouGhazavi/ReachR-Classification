clc
clear all

% Load the workspace 

% Dataset: R22-20_analysis_2023

% Data = load('/Users/Ines/Dropbox/Retinotectal_analysis_for Nilou/matlab_analysis/R22-20_analysis_2023/C1_YASS/C1_workspace.mat');

% Dataset:R22-22_analysis_2023
% Data = load('/Users/Ines/Dropbox/Retinotectal_analysis_for Nilou/matlab_analysis/R22-22_analysis_2023/C1_YASS/workspace.mat');

% % Dataset:R22-27_analysis_2023
Data = load('/Users/Ines/Dropbox/Retinotectal_analysis_for Nilou/matlab_analysis/R22-27_analysis_2023/C1_YASS/workspace.mat');

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

% Initialize cell array to store results for each cell
cell_s = cell(num_cells, 1);

% Initialize arrays to store ReachR and non-ReachR cells
ReachR_cells = [];
Non_ReachR_cells = [];


% Loop over each cell 
% 
for i=1:num_cells
    
    % spikes of each cell 
    cell_spikes = [];
    spike_times_epoch=[];
    % compare dimmest and brightest epochs 
    dimmest_epoch=[];
    brightest_epoch=[];

    % all the spikes for each cell 
    cell_spikes = spike_times{i};
    spike_times_epoch = cell_spikes;
    
  
     

    % find sampling frequency for each cell
    
    spike_times_cell1=spike_times{i,1};
    % calculate time intervals between consecutive spikes
    time_intervals = diff(spike_times_cell1);

   % calculate mean time interval (in seconds)
    mean_interval = mean(time_intervals);

   % calculate estimated sampling rate (in Hz)
   sampling_frequency=1 / mean_interval;
   cell_s{i}.fs = sampling_frequency;
    
   
   
    % trigger times for each epoch 
    dimmest_epoch=epoch_times{1};
    brightest_epoch=epoch_times{9};
    

    % Loop over each trigger in each epoch 
    spike_ratio_D = [];
    spike_ratio_B = [];
        
        


    % for the dimmest ( epoch=1)
    dimmest_spikes = {};
      
    
        % evoked response (response-baseline)
        evoked_response_D_vec{i}= 0;
        
        
        
        % epoch 8 or 1 
        for j = 1:length(epoch_times{1})

            spike_times_before_flash_D =[];
            spike_times_after_flash_D=[];

            spike_counts_before_D = [];
            spike_counts_after_D =[];
            evoked_response_D=[];
            
            % Each trigger in epoch=1, or 8 
            trig_time_epoch_D = epoch_times{1}(j);
            
            % Choose a narrow window (s) 5 ms or 10 ms ( after the trigger
            window=0.2;
            % Choose a window for baseline ( after the trigger)
            window_b=0.2;
            
            spike_times_before_flash_D = spike_times_epoch(spike_times_epoch >= trig_time_epoch_D - window_b & spike_times_epoch < trig_time_epoch_D);
            spike_times_after_flash_D = spike_times_epoch(spike_times_epoch > trig_time_epoch_D & spike_times_epoch < trig_time_epoch_D + window);
            evoked_response_D=length(spike_times_epoch(spike_times_epoch > trig_time_epoch_D & spike_times_epoch < trig_time_epoch_D + window)) - length(spike_times_epoch(spike_times_epoch >= trig_time_epoch_D - window_b & spike_times_epoch < trig_time_epoch_D));
            dimmest_spikes=spike_times_after_flash_D;
            
            % save them in cell_s struct ( apike activity after flash )
            cell_s{i}.cell_ID=cell_ids(i);
            cell_s{i}.epoch(1).triggered(j).spike_activity_before_D=spike_times_before_flash_D;
            cell_s{i}.epoch(1).triggered(j).spike_activity_after_D=spike_times_after_flash_D;
            cell_s{i}.epoch(1).triggered(j).evoked_response_D=evoked_response_D;
           
            
            
            % Add the current evoked response to the total 
            evoked_response_D_vec{i} =evoked_response_D+evoked_response_D_vec{i};
            
            
            % Method 1: Ratio of spike rate/spike counts before and after flash 

            spike_counts_before_D = length(spike_times_before_flash_D);
            spike_counts_after_D = length(spike_times_after_flash_D);

            % save them in cell_s struct ( apike activity after flash )
            
            cell_s{i}.epoch(1).triggered(j).spike_count_before_D= spike_counts_before_D;
            cell_s{i}.epoch(1).triggered(j).spike_count_after_D= spike_counts_after_D;
            
            % if there was no spike either before or after the flash 
            if spike_counts_before_D==0||spike_counts_after_D==0 

                if spike_counts_after_D==0 && spike_counts_before_D==0
                    
                    cell_s{i}.epoch(1).triggered(j).spike_ratio_D =1;

                elseif spike_counts_before_D==0 

                    cell_s{i}.epoch(1).triggered(j).spike_ratio_D = spike_counts_after_D ./ 1;
                else
                    
                    cell_s{i}.epoch(1).triggered(j).spike_ratio_D = 1./ spike_counts_before_D;
                end 
                 
            else 

                cell_s{i}.epoch(1).triggered(j).spike_ratio_D = spike_counts_after_D ./ spike_counts_before_D;
            
            end
        
        end
      
    
 
        
        
        % save the mean of the evoked response 
        cell_s{i}.epoch(1).mean_evoked_response_D=mean([cell_s{i}.epoch(1).triggered.evoked_response_D]);

        % save it 
        avg_D_cell(i)=cell_s{i}.epoch(1).mean_evoked_response_D;

        % find the average ( divide it by number of triggers in epoch 1) 
        avg_D=(evoked_response_D_vec{i})/length(epoch_times{1});   
    
        
        
        
    % For the brightest light (epoch=9)
        brightest_spikes={};
        evoked_response_B_vec{i}= 0;

        for j = 1:length(epoch_times{9})



            spike_times_before_flash_B =[];
            spike_times_after_flash_B=[];

            spike_counts_before_B = [];
            spike_counts_after_B =[];
            

            % Each trigger for epoch=9
            trig_time_epoch_B = epoch_times{9}(j);

            % Choose a very narrow window 5 ms or 10 ms 
            spike_times_before_flash_B = spike_times_epoch(spike_times_epoch >= trig_time_epoch_B - window_b & spike_times_epoch < trig_time_epoch_B);
            spike_times_after_flash_B = spike_times_epoch(spike_times_epoch > trig_time_epoch_B & spike_times_epoch < trig_time_epoch_B + window);
            evoked_response_B=length(spike_times_epoch(spike_times_epoch > trig_time_epoch_B & spike_times_epoch < trig_time_epoch_B + window)) - length(spike_times_epoch(spike_times_epoch >= trig_time_epoch_B - window_b & spike_times_epoch < trig_time_epoch_B));

            brightest_spikes=spike_times_after_flash_B;
            
            % save the spiking activity after flash 
            cell_s{i}.epoch(9).triggered(j).spike_activity_before_B=spike_times_before_flash_B;
            cell_s{i}.epoch(9).triggered(j).spike_activity_after_B=spike_times_after_flash_B;
            cell_s{i}.epoch(9).triggered(j).evoked_response_B=evoked_response_B;




            evoked_response_B_vec{i} =evoked_response_B+evoked_response_B_vec{i};

            % find number of spikes befoer and after flash 
            spike_counts_before_B = length(spike_times_before_flash_B);
            spike_counts_after_B = length(spike_times_after_flash_B);

            % save number of spikes after flash 
            
            
            cell_s{i}.epoch(9).triggered(j).spike_count_before_B=spike_counts_before_B;
            cell_s{i}.epoch(9).triggered(j).spike_count_after_B=spike_counts_after_B;
            
            evoked_response_B_vec{j} =evoked_response_B;
            
            if spike_counts_before_B==0|| spike_counts_after_B==0 

                if spike_counts_after_B==0 && spike_counts_before_B==0
                    
                    cell_s{i}.epoch(9).triggered(j).spike_ratio_B =1;

                elseif spike_counts_before_B==0 

                cell_s{i}.epoch(9).triggered(j).spike_ratio_B = spike_counts_after_B ./ 1;
                else
                    
                cell_s{i}.epoch(9).triggered(j).spike_ratio_B = 1./ spike_counts_before_B;

                end 
                 
            else 

                cell_s{i}.epoch(9).triggered(j).spike_ratio_B = spike_counts_after_B ./ spike_counts_before_B;
            
            end
            
            
        end


cell_s{i}.epoch(9).mean_evoked_response_B=mean([cell_s{i}.epoch(9).triggered.evoked_response_B]);
% using mean function 
avg_B_cell(i)=cell_s{i}.epoch(9).mean_evoked_response_B;

avg_B=(evoked_response_B_vec{i})/length(epoch_times{9}); 
    
    


end
save('spike_average_epoch9.mat','avg_B_cell')
save('spike_average_epoch1.mat','avg_D_cell')





%% Classify Different groups: ReaChR, Non-ReaCHR, Ambiguous

mask1 = (avg_B_cell > 2) & (avg_D_cell >= min(avg_D_cell)) & (avg_D_cell <=max(avg_D_cell));
mask2 = (avg_B_cell > 1) & (avg_B_cell < 2)& (avg_D_cell >= min(avg_D_cell)) & (avg_D_cell <=max(avg_D_cell));
mask3 = (avg_B_cell < 1)& (avg_D_cell >= min(avg_D_cell)) & (avg_D_cell <=max(avg_D_cell));

% Cell numbers
ReaChR_algorithm_cell_num=cells(mask1)
Ambiguous_algorithm_cell_num=cells(mask2)
Non_ReaChR_algorithm_cell_num=cells(mask3)

% Cell IDs 
ReaChR_algorithm_cell_ID=cell_ids(ReaChR_algorithm_cell_num)
Ambiguous_algorithm_cell_ID=cell_ids(Ambiguous_algorithm_cell_num)
Non_ReaChR_algorithm_cel_ID=cell_ids(Non_ReaChR_algorithm_cell_num)



% get a 2D plot
figure()
scatter(avg_D_cell(mask1), avg_B_cell(mask1), 100, 'g', 'filled');
cell_numbers_str = arrayfun(@(x) num2str(ReaChR_algorithm_cell_num(x)), 1:numel(ReaChR_algorithm_cell_num), 'UniformOutput', false);
text(avg_D_cell(ReaChR_algorithm_cell_num), avg_B_cell(ReaChR_algorithm_cell_num), cell_numbers_str, 'Color', 'black', 'FontSize', 10);

hold on
scatter(avg_D_cell(mask2), avg_B_cell(mask2), 100, 'r', 'filled');
cell_numbers_str = arrayfun(@(x) num2str(Ambiguous_algorithm_cell_num(x)), 1:numel(Ambiguous_algorithm_cell_num), 'UniformOutput', false);
text(avg_D_cell(Ambiguous_algorithm_cell_num), avg_B_cell(Ambiguous_algorithm_cell_num), cell_numbers_str, 'Color', 'black', 'FontSize', 10);
hold on 
scatter(avg_D_cell(mask3), avg_B_cell(mask3),100, 'b', 'filled');
cell_numbers_str = arrayfun(@(x) num2str(Non_ReaChR_algorithm_cell_num(x)), 1:numel(Non_ReaChR_algorithm_cell_num), 'UniformOutput', false);
text(avg_D_cell(Non_ReaChR_algorithm_cell_num), avg_B_cell(Non_ReaChR_algorithm_cell_num), cell_numbers_str, 'Color', 'black', 'FontSize', 10);

xlabel('Epoch 1 ');
ylabel('Epoch 9 (brightest)');
zlabel('Probability');
legend('ReaChR','Ambiguous','Non-ReaChR')









%% identify ReaCh R cells- Visualize the cells ( Comparison with visually identified cells)

% Load the dataset [ Remove this ] 

% Visually Identified cells 

% Visually inspected Dataset #1
ReaChR_IDv_1=load('/Users/Ines/Dropbox/Retinotectal_analysis_for Nilou/matlab_analysis/R22-20_analysis_2023/C1_YASS/ReaChR_ID.mat');
ReaChR_Chv_1=load('/Users/Ines/Dropbox/Retinotectal_analysis_for Nilou/matlab_analysis/R22-20_analysis_2023/C1_YASS/ReaChR_IDX.mat');

% Visually inspected Dataset #2
ReaChR_IDv_2=load('/Users/Ines/Dropbox/Retinotectal_analysis_for Nilou/matlab_analysis/R22-22_analysis_2023/C1_YASS/ReaChR_ID.mat');
ReaChR_Chv_2=load('/Users/Ines/Dropbox/Retinotectal_analysis_for Nilou/matlab_analysis/R22-22_analysis_2023/C1_YASS/ReaChR_IDX.mat');

% Visually inspected Dataset #3
ReaChR_IDv_3=load('/Users/Ines/Dropbox/Retinotectal_analysis_for Nilou/matlab_analysis/R22-27_analysis_2023/C1_YASS/ReaChR_ID.mat');
ReaChR_Chv_3=load('/Users/Ines/Dropbox/Retinotectal_analysis_for Nilou/matlab_analysis/R22-27_analysis_2023/C1_YASS/ReaChR_IDX.mat');


% Store the cell IDs in a vector for visually identified cells 
visually_identified_reachR_cells_1 = (ReaChR_IDv_1.reachr_id);
visually_identified_reachR_cells_2=(ReaChR_IDv_2.reachr_id);
visually_identified_reachR_cells_3=(ReaChR_IDv_3.reachr_id);

% cell IDs
Visually_identified=visually_identified_reachR_cells_3;
% convert cell IDs to cell numbers 
[~, idx] = (ismember(Visually_identified, cell_ids))



% Algorithm identified cells 

% Identify ReaChR cells 

x_edges = -20:0.1:20;
y_edges = -20:0.1:20;
bin_edges = {x_edges, y_edges};
cells=1:length(cell_ids);

% Compute histograms for each epoch
hist_D = histcounts2(cells, avg_D_cell, x_edges, y_edges, 'Normalization', 'probability');
hist_B = histcounts2(cells, avg_B_cell, x_edges, y_edges, 'Normalization', 'probability');

% Create a mask to select cells with epoch 9 > 5 and epoch 1 within -5 and
% 0 ( specify the boundary)
mask = (avg_B_cell > 2) & (avg_D_cell >= min(avg_D_cell)) & (avg_D_cell <=max(avg_D_cell));
% mask = (avg_B_cell >2) & (avg_D_cell >= 2);


%Not masked ( Non-ReaChR cells)
mask_not = ~mask;
cell_numbers=cells(mask)
cell_numbers_non_reachR=cells(mask_not)
% Cell IDs for ReaCh_R cells 
Algorithm_identified_ReaChR=cell_ids(cell_numbers)
Algorithm_identified_NonReaChR=cell_ids(cell_numbers_non_reachR)

% Save ReaChR cells idenified by algorithm

save('ReachR_cells_computationally.mat', 'ReachR_cells');
save('Non_ReachR_cells_computationally.mat', 'Non_ReachR_cells');
% ReaChR_IDs_Computationally=load('ReachR_cells_computationally.mat');

% Common ReaChR cells 
common_reachR_cells=intersect(Visually_identified,Algorithm_identified_ReaChR);
save('common_reachR_cells.mat','common_reachR_cells')

True_positive=intersect(Visually_identified,Algorithm_identified_ReaChR);
False_positive=setdiff(Visually_identified,Algorithm_identified_ReaChR);

num_ReaChR=length(Visually_identified);
num_false_positive=length(False_positive);
num_true_positive=length(True_positive);

% Visualize the results 

% Plot the 3D histogram
[counts, ~] = hist3([avg_D_cell', avg_B_cell'], 'Edges', bin_edges, 'DataDensity', 'on');

figure;
hist3([avg_D_cell', avg_B_cell'],'Edges', bin_edges, 'FaceColor', 'interp', 'CDataMode', 'auto');
set(gcf, 'Renderer', 'OpenGL'); % for better performance
set(get(gca,'child'),'FaceAlpha',0.5,'EdgeAlpha',0.1); % adjust opacity
xlabel('Epoch 1 ');
ylabel('Epoch 9 (brightest)');
zlabel('Counts');
legend('data','Misclassified ReaChR','Non-ReaChR','ReaChR')

% hide the first entry in the legend
colormap(jet);
caxis([0 200]);
colorbar;


% Use mask 
hold on 
scatter3(avg_D_cell(mask), avg_B_cell(mask), zeros(sum(mask), 1),50, 'r', 'filled');
hold on 
scatter3(avg_D_cell(mask_not), avg_B_cell(mask_not), zeros(sum(mask_not), 1), 50, 'b', 'filled');
hold on 
scatter3(avg_D_cell(idx), avg_B_cell(idx), zeros(length(idx),1), 50, 'g', 'filled');
legend('data','Misclassified ReaChR','Non-ReaChR','ReaChR')
hold on 

% Set the view to top-down
view([0 90]);

% get a 2D plot
figure()
scatter(avg_D_cell(mask), avg_B_cell(mask), 100, 'r', 'filled');
cell_numbers_str = arrayfun(@(x) num2str(cell_numbers(x)), 1:numel(cell_numbers), 'UniformOutput', false);
text(avg_D_cell(cell_numbers), avg_B_cell(cell_numbers), cell_numbers_str, 'Color', 'black', 'FontSize', 10);

hold on
scatter(avg_D_cell(mask_not), avg_B_cell(mask_not), 100, 'b', 'filled');
hold on 
scatter(avg_D_cell(idx), avg_B_cell(idx),100, 'g', 'filled');
common_cells = intersect(visually_identified_reachR_cells_2, cell_ids);
[~, common_indices] = ismember(common_cells, cell_ids);
cell_numbers = find(mask);
cell_number_str = arrayfun(@(x) num2str(idx(x)), 1:numel(idx), 'UniformOutput', false);
text(avg_D_cell(idx), avg_B_cell(idx), cell_number_str, 'Color', 'black', 'FontSize', 10);

xlabel('Epoch 1 ');
ylabel('Epoch 9 (brightest)');
zlabel('Probability');
legend('Misclassified ReaChR','Non-ReaChR','ReaChR')



% % Save the results and compare the results of the algorithm with visually identified ReaChR+

save('ReachR_cells_S.mat', 'ReachR_cells');
save('Non_ReachR_cells_S.mat', 'Non_ReachR_cells');







%% Load KNN model and make prediction for the random dataset 
loadedModel = load('knnModel.mat');
knnModel = loadedModel.knnModel;

% update the data 
newData=[avg_D_cell', avg_B_cell']

% predict the labels 
predictedLabelsNewData = predict(knnModel, newData);
% Plot the data points
figure();
gscatter(newData(:, 1), newData(:, 2), predictedLabelsNewData, 'rb', 'o', 10);

% Customize the plot
xlabel('Feature 1');
ylabel('Feature 2');
title('KNN Predictions on New Dataset');
legend('ReaChR', 'Non-ReaChR', 'Location', 'best');




%% Save the predicted labels 

num_rows=length(newData);
KNN_predicted_labels=cell(num_rows,4);
column_labels = {'Cell Number', 'Cell ID', 'Predicted Label', 'Predicted Label-Binary'};
KNN_predicted_labels = [column_labels; KNN_predicted_labels];
% Fill the KNN_predicted_labels cell array
for i = 1:num_rows
    KNN_predicted_labels{i, 1} = cells(i); % Assign the cell number
    KNN_predicted_labels{i, 2} = cell_ids(i); % Assign the cell IDs
    
    if predictedLabelsNewData(i) == 1
        KNN_predicted_labels{i, 3} = 'ReaChR';
        KNN_predicted_labels{i, 4} = 1;
        
    else
        KNN_predicted_labels{i, 3} = 'Non_ReaChR';
        KNN_predicted_labels{i, 4} = 0;
    end
end

ReaChR_KNN_cell_numbers = cells(predictedLabelsNewData == 1);
ReaChR_KNN_cell_IDs = cell_ids(predictedLabelsNewData == 1);

NonReaChR_KNN_cell_numbers = cells(predictedLabelsNewData == 0);
NonReaChR_KNN_cell_IDs = cell_ids(predictedLabelsNewData == 0);





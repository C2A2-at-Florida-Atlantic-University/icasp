% This script calculates inter-antenna offsets across frames of a 
% reference antenna across rows of the MIMO receiver.
%
% The output of this script will be stored in the upper-level folder
% 'data'. The file name is 'iao.mat'.

clear; close all;

% Immutable configurations
ANT_DIST_ROW = 0.07935;
ANT_COUNT = 4;
FREQ = 3.55e9;
FRAME_LEN = 2048; % defines how many IQ samples are there in each frame
LTS_LEN = 64; % defines how many IQ samples are there in each LTS symbol
ANT_ROWS = [9, 11, 13, 15; 17, 19, 21, 23; 1, 3, 5, 7; 41, 43, 45, 47; 25, 27, 29, 31; 33, 35, 37, 39]; % 2D array defining physical indexing of antennas

% Dataset configurations
MIMO_DATASETS = ["run1x1", "run2x1", "run3x1", "run4x1", "run5x1", "run6x1", "run7x1", "run8x1"];
FRAME_COUNT = 20;
LTS_COUNT = 1;
PATH_SOURCE = '../../../data/';
PATH_TARGET = '../../../data/iao.mat';

% Export configurations

% Modifiable configurations
LTS_START_SHIFT = 0; % defines the start of a "window" of LTS symbols that's used for estimation. is a multiplier of LTS_LEN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Go through each available dataset

% 1.1. Prepare a matrix to store inter-antenna offsets for all frames
frame_antenna_offsets = zeros(length(MIMO_DATASETS), FRAME_COUNT, 6, 3);

for file_idx = 1:length(MIMO_DATASETS)
   
    % 2. Load a single dataset
    mimo_run_data = load(strcat(PATH_SOURCE, MIMO_DATASETS(file_idx), ".mat"));
    
    % 2.1. Calculate ground truth angles
    [azimuth, elevation] = utils.get_ground_truth(getfield(mimo_run_data, 'client_coord'));
    
    % 3. Go through each of the required frames
    for frame_idx = 1:FRAME_COUNT
        frame_rmses = [];
    
        % 3.2. Go through each of the antenna rows, record inter-antenna offsets
        for row_idx = 1:6
            try
                % Prepare frame & antenna indexes we'll be working with
                frame = getfield(mimo_run_data, strcat('reference_frame_', int2str(frame_idx + 100)));
                reference_lts_positions = getfield(mimo_run_data, 'reference_lts_positions');
            
                % Prepare indexes of antennas in 1st row (w/ 1-based indexing)
                ant_row_idx = ANT_ROWS(row_idx, :) + 1;
            
                % Prepare indexes of LTS1 for each antenna in the row (w/ 1-based indexing)
                lts_indexes = reference_lts_positions(frame_idx + 1, ant_row_idx) + 1;
    
                % Add an offset (multiplier of LTS_LEN) to start indexes to implement a shifting window
                lts_indexes = lts_indexes + LTS_LEN * LTS_START_SHIFT;
            
                % Extract samples for each of the antennas we need
                ant1 = frame(ant_row_idx(1), :)';
                ant2 = frame(ant_row_idx(2), :)';
                ant3 = frame(ant_row_idx(3), :)';
                ant4 = frame(ant_row_idx(4), :)';
            
                % Correct phases within each of antennas using LTS1 & LTS2
                ant1_corr = utils.run_doppler_calibration(ant1, FRAME_LEN, lts_indexes(1), 'Antenna 1', LTS_LEN, false);
                ant2_corr = utils.run_doppler_calibration(ant2, FRAME_LEN, lts_indexes(2), 'Antenna 2', LTS_LEN, false);
                ant3_corr = utils.run_doppler_calibration(ant3, FRAME_LEN, lts_indexes(3), 'Antenna 3', LTS_LEN, false);
                ant4_corr = utils.run_doppler_calibration(ant4, FRAME_LEN, lts_indexes(4), 'Antenna 4', LTS_LEN, false);
            
                % Correct phases across antennas using lts_count LTS symbols
                [ant1_final, ant2_final, ant3_final, ant4_final, w_ant12, w_ant13, w_ant14] = utils.run_interant_calibration4( ...
                    lts_indexes, ant1_corr, ant2_corr, ant3_corr, ant4_corr, LTS_LEN, LTS_COUNT);
                
                % Record inter-antenna offsets
                frame_antenna_offsets(file_idx, frame_idx, row_idx, 1) = w_ant12;
                frame_antenna_offsets(file_idx, frame_idx, row_idx, 2) = w_ant13;
                frame_antenna_offsets(file_idx, frame_idx, row_idx, 3) = w_ant14;
            catch error
                frame_antenna_offsets(file_idx, frame_idx, row_idx, 1) = -1;
                frame_antenna_offsets(file_idx, frame_idx, row_idx, 2) = -1;
                frame_antenna_offsets(file_idx, frame_idx, row_idx, 3) = -1;
            end
        end
    end
end

% 4. Evaluate inter-antenna offsets

% 4.1. Prepare matrix for storing averaged (or median) offsets for further usage
median_offsets_final = zeros(length(MIMO_DATASETS), 6, 3);

for row_idx = 1:6
    f = figure('Name', strcat('Inter-antenna Phase Offsets. LTS: ', num2str(LTS_COUNT), '. ROW: ', num2str(row_idx)), 'NumberTitle', 'off');
    f.Position = [100 100 2000 800];
    for i = 1:size(frame_antenna_offsets, 1)
        % Calculate average offsets for each antenna pair
        offset_12 = median(unwrap(frame_antenna_offsets(i, :, row_idx, 1)));
        offset_13 = median(unwrap(frame_antenna_offsets(i, :, row_idx, 2)));
        offset_14 = median(unwrap(frame_antenna_offsets(i, :, row_idx, 3)));
    
        % Plot everything
        subplot(2, 4, i)
        plot(1:FRAME_COUNT, unwrap(frame_antenna_offsets(i, :, row_idx, 1)), '-red', 'DisplayName', 'Antenna 1-2');
        hold on
        plot(1:FRAME_COUNT, unwrap(frame_antenna_offsets(i, :, row_idx, 2)), '-green', 'DisplayName', 'Antenna 1-3');
        plot(1:FRAME_COUNT, unwrap(frame_antenna_offsets(i, :, row_idx, 3)), '-blue','DisplayName', 'Antenna 1-4');
        plot([1, FRAME_COUNT-1], [offset_12, offset_12], 'DisplayName', 'Antenna 1-2 (median)')
        plot([1, FRAME_COUNT-1], [offset_13, offset_13], 'DisplayName', 'Antenna 1-3 (median)')
        plot([1, FRAME_COUNT-1], [offset_14, offset_14], 'DisplayName', 'Antenna 1-4 (median)')
        hold off
        legend
        grid on
        xlabel('Frame Index');
        ylabel('Phase Offset, rad');
        ylim([-3*pi, 2*pi]);
        title(strcat(MIMO_DATASETS(i), ' (unwrapped)'));
    
        % Save median offsets
        median_offsets_final(i, row_idx, 1) = offset_12;
        median_offsets_final(i, row_idx, 2) = offset_13;
        median_offsets_final(i, row_idx, 3) = offset_14;
    end
end

% Export median offset matrix
save(PATH_TARGET, 'median_offsets_final');
clear; close all;

% This script generates a final batch of POWDER measurements, ready for
% evaluating MUSIC and Hankel SVD. These measurements are calibrated on
% a row-by-row basis.
% 
% The files will be saved in the 'data' folder.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('../../../data/iao_rows.mat');

FRAME_LEN = 2048; % defines how many IQ samples are there in each frame
LTS_LEN = 64; % defines how many IQ samples are there in each LTS symbol
ANT_ROWS = [9, 11, 13, 15; 17, 19, 21, 23; 1, 3, 5, 7; 41, 43, 45, 47; 25, 27, 29, 31; 33, 35, 37, 39]; % 2D array defining physical indexing of antennas
MIMO_DATASETS = ["run1x1", "run2x1", "run3x1", "run4x1", "run5x1", "run6x1", "run7x1", "run8x1"];
FRAME_COUNT = 19;
LTS_COUNT = 2;

PATH_SOURCE = '../../../data/';
PATH_TARGET = '../../../data/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(MIMO_DATASETS)
    disp(strcat('Dataset ', int2str(i)));
    generate_data_reference(MIMO_DATASETS(i), i, PATH_SOURCE, PATH_TARGET, FRAME_COUNT, LTS_COUNT, FRAME_LEN, LTS_LEN, ANT_ROWS);
    generate_data_client(MIMO_DATASETS(i), i, PATH_SOURCE, PATH_TARGET, FRAME_COUNT, LTS_COUNT, FRAME_LEN, LTS_LEN, ANT_ROWS, median_offsets_final);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = generate_data_client(dataset_name, client_idx, path_source, path_target, frame_count, lts_count, frame_len, lts_len, ant_rows, iao)
    % 1. Load original dataset
    mimo_run_data = load(strcat(path_source, dataset_name, '.mat'));
    [azimuth, ~] = utils.get_ground_truth(getfield(mimo_run_data, 'client_coord'));
    
    % 2. Perform processing for each of the available 20 frames
    output_samples = zeros(6, 4, 128, 20);
    
    for frame_idx = 1:frame_count
        for row_idx = 1:6
            try
                % Prepare indexes of antennas in 1st row (w/ 1-based indexing)
                ant_row_idx = ant_rows(row_idx, :) + 1;
            
                % Retrieve frame samples
                client_frame = getfield(mimo_run_data, strcat('client_frame_', int2str(frame_idx + 100)));
    
                % Determine indexes where first LTS'es begin
                client_lts_positions = getfield(mimo_run_data, 'client_lts_positions');
                client_lts_indexes = client_lts_positions(frame_idx + 1, ant_row_idx) + 1;
        
                % Extract samples for each of the antennas we need
                client_ant1 = client_frame(ant_row_idx(1), :)';
                client_ant2 = client_frame(ant_row_idx(2), :)';
                client_ant3 = client_frame(ant_row_idx(3), :)';
                client_ant4 = client_frame(ant_row_idx(4), :)';
            
                % Correct phases within each of antennas using LTS1 & LTS2
                client_ant1_corr = utils.run_doppler_calibration(client_ant1, frame_len, client_lts_indexes(1), 'Antenna 1', lts_len, false);
                client_ant2_corr = utils.run_doppler_calibration(client_ant2, frame_len, client_lts_indexes(2), 'Antenna 2', lts_len, false);
                client_ant3_corr = utils.run_doppler_calibration(client_ant3, frame_len, client_lts_indexes(3), 'Antenna 3', lts_len, false);
                client_ant4_corr = utils.run_doppler_calibration(client_ant4, frame_len, client_lts_indexes(4), 'Antenna 4', lts_len, false);
        
                % Use inter-antenna offsets for reference emitter to correct client samples
                % Note: uses median inter-antenna phase offsets. Calculated as a
                % median of phase offsets across mutiple consecutive reference frames.
                client_ant1_final = client_ant1_corr;
                client_ant2_final = client_ant2_corr.*exp(-1j * iao(client_idx, row_idx, 1));
                client_ant3_final = client_ant3_corr.*exp(-1j * iao(client_idx, row_idx, 2));
                client_ant4_final = client_ant4_corr.*exp(-1j * iao(client_idx, row_idx, 3));
        
                % Prepare MUSIC input array (1D MUSIC can be run on this matrix)
                client_X = zeros(4, lts_len * lts_count); 
                client_X(1, :) = client_ant1_final(client_lts_indexes(1) : client_lts_indexes(1) + lts_len * lts_count - 1)';
                client_X(2, :) = client_ant2_final(client_lts_indexes(2) : client_lts_indexes(2) + lts_len * lts_count - 1)';
                client_X(3, :) = client_ant3_final(client_lts_indexes(3) : client_lts_indexes(3) + lts_len * lts_count - 1)';
                client_X(4, :) = client_ant4_final(client_lts_indexes(4) : client_lts_indexes(4) + lts_len * lts_count - 1)';
    
                % Prepare output (take only 1st sample from the array)
                output_samples(row_idx, 1, :, frame_idx) = reshape(client_X(1, :), [128, 1]);
                output_samples(row_idx, 2, :, frame_idx) = reshape(client_X(2, :), [128, 1]);
                output_samples(row_idx, 3, :, frame_idx) = reshape(client_X(3, :), [128, 1]);
                output_samples(row_idx, 4, :, frame_idx) = reshape(client_X(4, :), [128, 1]);
            catch error
                % disp(strcat('Problems with LTS indexing. Row #', int2str(row_idx)));
    
                % Prepare output (take only 1st sample from the array)
                output_samples(row_idx, 1, :, frame_idx) = NaN;
                output_samples(row_idx, 2, :, frame_idx) = NaN;
                output_samples(row_idx, 3, :, frame_idx) = NaN;
                output_samples(row_idx, 4, :, frame_idx) = NaN;
            end
        end
    end
    
    % 3. Save output array for further processing
    save_filename = strcat(path_target, 'client_x', num2str(client_idx), '_for_az.mat');
    for i = 1:frame_count
        varname = strcat('output_samples_frame_', num2str(i));
        S.(varname) = squeeze(output_samples(:, :, :, i));
    end
    S.output_az = azimuth;
    save(save_filename, '-struct', 'S');
end

function [] = generate_data_reference(dataset_name, client_idx, path_source, path_target, frame_count, lts_count, frame_len, lts_len, ant_rows)

    % 1. Load original dataset
    mimo_run_data = load(strcat(path_source, dataset_name, ".mat"));
    [azimuth, ~] = utils.get_ground_truth(getfield(mimo_run_data, 'ref_coord'));
    
    % 2. Perform processing for each of the available 20 frames
    output_samples = zeros(6, 4, 128, 20);
    
    for frame_idx = 1:frame_count
        for row_idx = 1:6
            try
                % Prepare frame & antenna indexes we'll be working with
                frame = getfield(mimo_run_data, strcat('reference_frame_', int2str(frame_idx + 100)));
                reference_lts_positions = getfield(mimo_run_data, 'reference_lts_positions');
            
                % Prepare indexes of antennas in 1st row (w/ 1-based indexing)
                ant_row_idx = ant_rows(row_idx, :) + 1;
            
                % Prepare indexes of LTS1 for each antenna in the row (w/ 1-based indexing)
                lts_indexes = reference_lts_positions(frame_idx + 1, ant_row_idx) + 1;
            
                % Extract samples for each of the antennas we need
                ant1 = frame(ant_row_idx(1), :)';
                ant2 = frame(ant_row_idx(2), :)';
                ant3 = frame(ant_row_idx(3), :)';
                ant4 = frame(ant_row_idx(4), :)';
            
                % Correct phases within each of antennas using LTS1 & LTS2
                ant1_corr = utils.run_doppler_calibration(ant1, frame_len, lts_indexes(1), 'Antenna 1', lts_len, false);
                ant2_corr = utils.run_doppler_calibration(ant2, frame_len, lts_indexes(2), 'Antenna 2', lts_len, false);
                ant3_corr = utils.run_doppler_calibration(ant3, frame_len, lts_indexes(3), 'Antenna 3', lts_len, false);
                ant4_corr = utils.run_doppler_calibration(ant4, frame_len, lts_indexes(4), 'Antenna 4', lts_len, false);
            
                % Correct phases across antennas using lts_count LTS symbols
                [ant1_final, ant2_final, ant3_final, ant4_final, ~, ~, ~] = utils.run_interant_calibration4( ...
                    lts_indexes, ant1_corr, ant2_corr, ant3_corr, ant4_corr, lts_len, lts_count);
    
                % Run MUSIC on these four antennas
                X = zeros(4, lts_len * lts_count); 
    
                % Try final calibrated data
                X(1, :) = ant1_final(lts_indexes(1) : lts_indexes(1) + lts_len * lts_count - 1)';
                X(2, :) = ant2_final(lts_indexes(2) : lts_indexes(2) + lts_len * lts_count - 1)';
                X(3, :) = ant3_final(lts_indexes(3) : lts_indexes(3) + lts_len * lts_count - 1)';
                X(4, :) = ant4_final(lts_indexes(4) : lts_indexes(4) + lts_len * lts_count - 1)';
    
                % Prepare output (take only 1st sample from the array)
                output_samples(row_idx, 1, :, frame_idx) = reshape(X(1, :), [128, 1]);
                output_samples(row_idx, 2, :, frame_idx) = reshape(X(2, :), [128, 1]);
                output_samples(row_idx, 3, :, frame_idx) = reshape(X(3, :), [128, 1]);
                output_samples(row_idx, 4, :, frame_idx) = reshape(X(4, :), [128, 1]);
            catch error
                % disp(strcat('Problems with LTS indexing. Row #', int2str(row_idx)));
    
                % Prepare output (take only 1st sample from the array)
                output_samples(row_idx, 1, :, frame_idx) = NaN;
                output_samples(row_idx, 2, :, frame_idx) = NaN;
                output_samples(row_idx, 3, :, frame_idx) = NaN;
                output_samples(row_idx, 4, :, frame_idx) = NaN;
            end
        end
    end
    
    % 3. Save output array for further processing
    save_filename = strcat(path_target, 'reference_x', num2str(client_idx), '_for_az.mat');
    for i = 1:frame_count
        varname = strcat('output_samples_frame_', num2str(i));
        S.(varname) = squeeze(output_samples(:, :, :, i));
    end
    
    S.output_az = azimuth;
    save(save_filename, '-struct', 'S');
end
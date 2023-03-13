clear all; clc; close all; 

addpath('../utilities/');

PATH_SOURCE = '../../../../data/';

rand_seed = 100;

doa_true_angles_az = [-15.44, -6.93, 4.84, 17.72, 27.88, 36.93, 53, 76, 0];
doa_true_angles_elev = [-50, -49.1, -48.1, -48.1, -50, -43.07, -34.16, -29.25, -85.67];
T = 128; % number of samples to simulate
f = 3.55e9; % frequency
lambda = 0.0844486; % wavelength (meters)
d_row = 0.07935; % antenna spacing within a row (meters)
d_col = 0.06668; % antenna spacing within a column (meters)
SNR = 1; % SNR(dB)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doa_threshold = 10;

% full_test_sim_elev(doa_true_angles_elev, d_col, T, lambda, 6, 20, SNR, rand_seed, doa_threshold);
full_test_sim_az(doa_true_angles_az, d_row, T, lambda, 4, 20, SNR, rand_seed, doa_threshold);

% full_test_real_elev(PATH_SOURCE, lambda, 6, 19, doa_threshold);
% full_test_real_az(PATH_SOURCE, lambda, 4, 19, doa_threshold);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = full_test_real_elev(path_source, lambda, Nr, N, doa_threshold)
    % 1. Read measurements
    [doa_true_angles_elev, real_signal_elev] = read_real_data_elev(path_source);
    
    % 2. Filter samples in each client & frame which are below given DoA threshold
    signal_real_elev_filtered = filter_signal(real_signal_elev, doa_true_angles_elev, Nr, doa_threshold);
    
    % 3. Run Hankel vs. MUSIC evaluation on filtered dataset
    rmse_hankel_elev_filtered = evaluate_hankel(doa_true_angles_elev, signal_real_elev_filtered, Nr, N);
    rmse_music_elev_filtered = evaluate_music(doa_true_angles_elev, signal_real_elev_filtered, lambda, Nr);

    % 4. Run Hankel vs. MUSIC evaluation on full dataset (as comparison)
    signal_real_elev_full = convert_signal_to_cell(real_signal_elev);
    rmse_hankel_elev_full = evaluate_hankel(doa_true_angles_elev, signal_real_elev_full, Nr, N);
    rmse_music_elev_full = evaluate_music(doa_true_angles_elev, signal_real_elev_full, lambda, Nr);
    
    % Plot Hankel performance
    figure
    hold on;
    plot(1:9, rmse_hankel_elev_filtered, 'red-o','MarkerFaceColor', 'red', 'DisplayName','Hankel SVD Filtered');
    plot(1:9, rmse_music_elev_filtered, 'green-o','MarkerFaceColor', 'green', 'DisplayName','MUSIC Filtered');
    plot(1:9, rmse_hankel_elev_full, 'yellow-o','MarkerFaceColor', 'yellow', 'DisplayName','Hankel SVD Full');
    plot(1:9, rmse_music_elev_full, 'blue-o','MarkerFaceColor', 'blue', 'DisplayName','MUSIC Full');
    hold off;
    grid on; 
    xlabel('Client');
    ylabel('DOA estimate RMSE');
    legend();
    title("Hankel vs. MUSIC. Elevation. Real Data.");
end

function [] = full_test_real_az(path_source, lambda, Nr, N, doa_threshold)
    % 1. Read measurements
    [doa_true_angles_az, real_signal_az] = read_real_data_az(path_source);
    
    % 2. Filter samples in each client & frame which are below given DoA threshold
    test_signal_sim_az_filtered = filter_signal(real_signal_az, doa_true_angles_az, Nr, doa_threshold);
    
    % 3. Run Hanke vs. MUSIC evaluation on filtered dataset
    rmse_hankel_az_filtered = evaluate_hankel(doa_true_angles_az, test_signal_sim_az_filtered, Nr, N);
    rmse_music_az_filtered = evaluate_music(doa_true_angles_az, test_signal_sim_az_filtered, lambda, Nr);

    % 4. Run Hankel vs. MUSIC evaluation on full dataset (as comparison)
    signal_real_az_full = convert_signal_to_cell(real_signal_az);
    rmse_hankel_az_full = evaluate_hankel(doa_true_angles_az, signal_real_az_full, Nr, N);
    rmse_music_az_full = evaluate_music(doa_true_angles_az, signal_real_az_full, lambda, Nr);
    
    % Plot Hankel performance
    figure
    hold on;
    plot(1:9, rmse_hankel_az_filtered, 'red-o','MarkerFaceColor', 'red', 'DisplayName','Hankel SVD Filtered');
    plot(1:9, rmse_music_az_filtered, 'green-o','MarkerFaceColor', 'green', 'DisplayName','MUSIC Filtered');
    plot(1:9, rmse_hankel_az_full, 'yellow-o','MarkerFaceColor', 'yellow', 'DisplayName','Hankel SVD Full');
    plot(1:9, rmse_music_az_full, 'blue-o','MarkerFaceColor', 'blue', 'DisplayName','MUSIC Full');
    hold off;
    grid on; 
    xlabel('Client');
    ylabel('DOA estimate RMSE');
    legend();
    title("Hankel vs. MUSIC. Azimuth. Real Data.");
end

function [] = full_test_sim_elev(doa_true_angles_elev, d_col, T, lambda, Nr, N, SNR, rand_seed, doa_threshold)
    % 1. Generate simulated measurements
    sim_signal_elev = simulate_data(doa_true_angles_elev, T, lambda, d_col, Nr, SNR, N, rand_seed);
    
    % 2. Filter samples in each client & frame which are below given DoA threshold
    % Create array of cells of a size (M x N), where M - # of clients, N - # of frames
    test_signal_sim_elev_filtered = filter_signal(sim_signal_elev, doa_true_angles_elev, Nr, doa_threshold);
    
    % 3. Run Hanke vs. MUSIC evaluation on filtered dataset
    rmse_hankel_elev_filtered = evaluate_hankel(doa_true_angles_elev, test_signal_sim_elev_filtered, Nr, N);
    rmse_music_elev_filtered = evaluate_music(doa_true_angles_elev, test_signal_sim_elev_filtered, lambda, Nr);

    % 4. Run Hankel vs. MUSIC evaluation on full dataset (as comparison)
    signal_sim_elev_full = convert_signal_to_cell(sim_signal_elev);
    rmse_hankel_elev_full = evaluate_hankel(doa_true_angles_elev, signal_sim_elev_full, Nr, N);
    rmse_music_elev_full = evaluate_music(doa_true_angles_elev, signal_sim_elev_full, lambda, Nr);
    
    % Plot Hankel performance
    figure
    hold on;
    plot(1:9, rmse_hankel_elev_filtered, 'red-o','MarkerFaceColor', 'red', 'DisplayName','Hankel SVD Filtered');
    plot(1:9, rmse_music_elev_filtered, 'green-o','MarkerFaceColor', 'green', 'DisplayName','MUSIC Filtered');
    plot(1:9, rmse_hankel_elev_full, 'yellow-o','MarkerFaceColor', 'yellow', 'DisplayName','Hankel SVD Full');
    plot(1:9, rmse_music_elev_full, 'blue-o','MarkerFaceColor', 'blue', 'DisplayName','MUSIC Full');
    hold off;
    grid on; 
    xlabel('Client');
    ylabel('DOA estimate RMSE');
    legend();
    title("Hankel vs. MUSIC. Elevation. Simulated Data.");
end

function [] = full_test_sim_az(doa_true_angles_az, d_row, T, lambda, Nr, N, SNR, rand_seed, doa_threshold)
    % 1. Generate simulated measurements
    sim_signal_az = simulate_data(doa_true_angles_az, T, lambda, d_row, Nr, SNR, N, rand_seed);
    
    % 2. Filter samples in each client & frame which are below given DoA threshold
    % Create array of cells of a size (M x N), where M - # of clients, N - # of frames
    test_signal_sim_az_filtered = filter_signal(sim_signal_az, doa_true_angles_az, Nr, doa_threshold);
    
    % 3. Run Hanke vs. MUSIC evaluation on filtered dataset
    rmse_hankel_az_filtered = evaluate_hankel(doa_true_angles_az, test_signal_sim_az_filtered, Nr, N);
    rmse_music_az_filtered = evaluate_music(doa_true_angles_az, test_signal_sim_az_filtered, lambda, Nr);

    % 4. Run Hankel vs. MUSIC evaluation on full dataset (as comparison)
    signal_sim_az_full = convert_signal_to_cell(sim_signal_az);
    rmse_hankel_az_full = evaluate_hankel(doa_true_angles_az, signal_sim_az_full, Nr, N);
    rmse_music_az_full = evaluate_music(doa_true_angles_az, signal_sim_az_full, lambda, Nr);
    
    % Plot Hankel performance
    figure
    hold on;
    plot(1:9, rmse_hankel_az_filtered, 'red-o','MarkerFaceColor', 'red', 'DisplayName','Hankel SVD Filtered');
    plot(1:9, rmse_music_az_filtered, 'green-o','MarkerFaceColor', 'green', 'DisplayName','MUSIC Filtered');
    plot(1:9, rmse_hankel_az_full, 'yellow-o','MarkerFaceColor', 'yellow', 'DisplayName','Hankel SVD Full');
    plot(1:9, rmse_music_az_full, 'blue-o','MarkerFaceColor', 'blue', 'DisplayName','MUSIC Full');
    hold off;
    grid on; 
    xlabel('Client');
    ylabel('DOA estimate RMSE');
    legend();
    title("Hankel vs. MUSIC. Azimuth. Simulated Data.");
end

function [signal_filtered] = filter_signal(signal, true_doa, Nr, doa_threshold)
    signal_filtered = cell(size(signal, 3), size(signal, 4));
    for client_idx = 1:size(signal_filtered, 1)
        fprintf(strcat('\n\nClient: ', int2str(client_idx)));
        for frame_idx = 1:size(signal_filtered, 2)
            fprintf(strcat('\nFrame: ', int2str(frame_idx)));
    
            % Determine which samples within this frame are "within threshold"
            client_true_doa = true_doa(client_idx);
            client_samples = squeeze(signal(:, :, client_idx, frame_idx));
    
            frame_good_samples = find_frame_good_samples_hankel(client_true_doa, client_samples, Nr, doa_threshold, client_idx, false);
    
            % Extract "good" samples from the dataset & save filtered samples into the cell array
            signal_filtered(client_idx, frame_idx) = {client_samples(:, frame_good_samples)};
        end
    end
end

function [signal_full] = convert_signal_to_cell(signal)
    signal_full = cell(size(signal, 3), size(signal, 4));
    for client_idx = 1:size(signal_full, 1)
        for frame_idx = 1:size(signal_full, 2)
            signal_full(client_idx, frame_idx) = {squeeze(signal(:, :, client_idx, frame_idx))};
        end
    end
end

% Determines which samples in a frame produce sufficient accuracy based on threshold
function [sample_good_idx] = find_frame_good_samples_hankel(true_doa, frame_samples, Nr, doa_threshold, client_idx, plot_errors)
    frame_errors = zeros(128, 1);
    for sample_idx = 1:size(frame_errors, 1)
        X = squeeze(frame_samples(:, sample_idx));

        [svd_final_doa, ~, ~, ~, ~] = MSSA_SVD_findOptimalWandK_filtering_onetime_max_magnitude(X, Nr, true_doa, 2, 3, 1, 's', 1);

        % [mf_doa,~, ~] = aoa(X, Nr, 1);

        frame_errors(sample_idx, 1) = abs(true_doa - svd_final_doa);
    end
    sample_good_idx = frame_errors < doa_threshold;

    if plot_errors
        figure
        scatter(1:sum(sample_good_idx), frame_errors(sample_good_idx));
        title(strcat('Client', int2str(client_idx)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%% READING / SIMULATING DATA %%%%%%%%%%%%%%%%%%%%%%%%

function [test_signal_sim] = simulate_data(doa_true_angles, T, lambda, d_x, Nr, SNR, N, rand_seed)
    rng(rand_seed);
    test_signal_sim = zeros(Nr, T, length(doa_true_angles), N);
    for doa_idx = 1:length(doa_true_angles)
        for frame_idx = 1:N
            doa_true = doa_true_angles(doa_idx);
            doa_signal_sim = simulation.simulate_samples(doa_true, T, lambda, d_x, Nr, SNR);
            test_signal_sim(:, :, doa_idx, frame_idx) = doa_signal_sim;
        end
    end
end

function [doa_true_angles_az, real_signal_az] = read_real_data_az(path_source)
    real_signal_client = zeros(4, 128, 8, 19);
    real_signal_reference = zeros(4, 128, 1, 19);
    true_doa_client = zeros(8, 1);
    true_doa_reference = zeros(1, 1);
    for client_idx = 1:size(real_signal_client, 3)
        file_name = strcat(path_source, 'client_x' + string(client_idx) + '_for_az.mat'); 

        % Retrieve true DoA
        true_doa_client(client_idx, 1) = double(load(file_name).('output_az'));

        % Retrieve frames
        for frame_idx = 1:size(real_signal_client, 4)
            frame_name = 'output_samples_frame_' + string(frame_idx); 
            frame_matrix = load(file_name).(frame_name); %data matrix is a tensor (rows, cols, sample)

            for sample_idx = 1:128
                data_matrix = frame_matrix(:, :, sample_idx); 

                NaN_matrix = isnan(data_matrix); 

                [row_with_NaN,~] = find(NaN_matrix == 1); % find the rows that contain NaNs 
                
                data_matrix(row_with_NaN,:) = []; % remove this row/rows from data_matrix 
  
                averaged_row = mean(data_matrix,1).';
                real_signal_client(:, sample_idx, client_idx, frame_idx) = averaged_row;
            end
        end
    end

    for reference_idx = 1:size(real_signal_reference, 3)
        file_name = strcat(path_source, 'reference_x' + string(reference_idx) + '_for_az.mat'); 

        % Retrieve true DoA
        true_doa_reference(reference_idx, 1) = double(load(file_name).('output_az'));

        % Retrieve frames
        for frame_idx = 1:size(real_signal_reference, 4)
            frame_name = 'output_samples_frame_' + string(frame_idx); 
            frame_matrix = load(file_name).(frame_name); %data matrix is a tensor (rows, cols, sample)

            for sample_idx = 1:128
                data_matrix = frame_matrix(:, :, sample_idx); 

                NaN_matrix = isnan(data_matrix); 

                [row_with_NaN,~] = find(NaN_matrix == 1); % find the rows that contain NaNs 
                
                data_matrix(row_with_NaN,:) = []; % remove this row/rows from data_matrix 
  
                averaged_row = mean(data_matrix,1).';
                real_signal_reference(:, sample_idx, reference_idx, frame_idx) = averaged_row;
            end
        end
    end

    doa_true_angles_az = zeros(9, 1);
    doa_true_angles_az(1:8, 1) = true_doa_client;
    doa_true_angles_az(9, 1) = true_doa_reference;

    real_signal_az = zeros(4, 128, 9, 19);
    real_signal_az(:, :, 1:8, :) = real_signal_client;
    real_signal_az(:, :, 9, :) = real_signal_reference;
end

function [doa_true_angles_elev, real_signal_elev] = read_real_data_elev(path_source)
    real_signal_client = zeros(6, 128, 8, 19);
    real_signal_reference = zeros(6, 128, 1, 19);
    true_doa_client = zeros(8, 1);
    true_doa_reference = zeros(1, 1);
    for client_idx = 1:size(real_signal_client, 3)
        file_name = strcat(path_source, 'client_x' + string(client_idx) + '_for_elev.mat'); 

        % Retrieve true DoA
        true_doa_client(client_idx, 1) = double(load(file_name).('output_elev'));

        % Retrieve frames
        for frame_idx = 1:size(real_signal_client, 4)
            frame_name = 'output_samples_frame_' + string(frame_idx); 
            frame_matrix = load(file_name).(frame_name); %data matrix is a tensor (rows, cols, sample)

            for sample_idx = 1:128
                data_matrix = frame_matrix(:, :, sample_idx); 

                NaN_matrix = isnan(data_matrix); 

                [~,col_with_NaN] = find(NaN_matrix == 1); % find the rows that contain NaNs 
                
                data_matrix(:,col_with_NaN) = []; % remove this row/rows from data_matrix 
  
                averaged_col = mean(data_matrix,2);
                real_signal_client(:, sample_idx, client_idx, frame_idx) = averaged_col;
            end
        end
    end

    for reference_idx = 1:size(real_signal_reference, 3)
        file_name = strcat(path_source, 'reference_x' + string(reference_idx) + '_for_elev.mat'); 

        % Retrieve true DoA
        true_doa_reference(reference_idx, 1) = double(load(file_name).('output_elev'));

        % Retrieve frames
        for frame_idx = 1:size(real_signal_reference, 4)
            frame_name = 'output_samples_frame_' + string(frame_idx); 
            frame_matrix = load(file_name).(frame_name); %data matrix is a tensor (rows, cols, sample)

            for sample_idx = 1:128
                data_matrix = frame_matrix(:, :, sample_idx); 

                NaN_matrix = isnan(data_matrix); 

                [~,col_with_NaN] = find(NaN_matrix == 1); % find the rows that contain NaNs 
                
                data_matrix(:,col_with_NaN) = []; % remove this row/rows from data_matrix 
  
                averaged_col = mean(data_matrix,2);
                real_signal_reference(:, sample_idx, reference_idx, frame_idx) = averaged_col;
            end
        end
    end

    doa_true_angles_elev = zeros(9, 1);
    doa_true_angles_elev(1:8, 1) = true_doa_client;
    doa_true_angles_elev(9, 1) = true_doa_reference;

    real_signal_elev = zeros(6, 128, 9, 19);
    real_signal_elev(:, :, 1:8, :) = real_signal_client;
    real_signal_elev(:, :, 9, :) = real_signal_reference;
end

%%%%%%%%%%%%%%%%%%%%%%%%% RUNNING THE METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [doa_rmse] = evaluate_music(doa_true_angles, test_signal_sim, lambda, Nr)
    doa_rmse = zeros(size(test_signal_sim, 1), 1);
    for doa_idx = 1:length(doa_true_angles)

        true_doa = doa_true_angles(doa_idx);

        frame_errors = zeros(size(test_signal_sim, 2), 1);
        for frame_idx = 1:size(test_signal_sim, 2)
            X = cell2mat(test_signal_sim(doa_idx, frame_idx));

            doa_estimated = music.run_music(X, 1, lambda, Nr);
            if isnan(doa_estimated)
                frame_errors(frame_idx, 1) = NaN;
            else
                frame_errors(frame_idx, 1) = (doa_estimated - true_doa)^2;
            end
        end
        doa_rmse(doa_idx, 1) = sqrt(mean(frame_errors, 'omitnan'));
    end
end

function [doa_rmse] = evaluate_hankel(doa_true_angles, test_signal_sim, Nr, N)
    doa_rmse = zeros(length(doa_true_angles), 1);
    
    for doa_idx = 1:length(doa_true_angles) 
        fprintf(strcat('\n\nClient', int2str(doa_idx)));
   
        true_doa = doa_true_angles(doa_idx);
        frame_rmse = zeros(N, 1);
        
        for frame_idx = 1:N
            fprintf(strcat('\nFrame', int2str(frame_idx), ':'));

            data_matrix = cell2mat(test_signal_sim(doa_idx, frame_idx));
            
            sample_rmse = zeros(size(data_matrix, 2), 1);
            for sample_idx = 1:size(data_matrix, 2)
                fprintf('.');
                X = data_matrix(:,sample_idx);
                [svd_final_aoa, ~, ~, ~, ~] = MSSA_SVD_findOptimalWandK_filtering_onetime_max_magnitude(X, Nr, true_doa, 2, 3, 1, 's', 1);
                sample_rmse(sample_idx, 1) = (svd_final_aoa - true_doa)^2;
            end 
            frame_rmse(frame_idx, 1) = sqrt(mean(sample_rmse, 'omitnan'));
        end
        doa_rmse(doa_idx, 1) = mean(frame_rmse, 'omitnan');
    end 
end

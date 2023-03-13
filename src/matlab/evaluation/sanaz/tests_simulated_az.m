clear all; clc; close all; 

addpath('../utilities/');

rand_seed = 100;

doa_true_angles_az = [-15.44, -6.93, 4.84, 17.72, 27.88, 36.93, 53, 76, 0];
doa_true_angles_elev = [-50, -49.1, -48.1, -48.1, -50, -43.07, -34.16, -29.25, -85.67];
T = 128; % number of samples to simulate
f = 3.55e9; % frequency
lambda = 0.0844486; % wavelength (meters)
d_row = 0.07935; % antenna spacing within a row (meters)
d_col = 0.06668; % antenna spacing within a column (meters)
Nr  = 4; % number of receiver's antennas 
SNR = 1; % SNR(dB)
N = 20; % number of frames to simulate for each test case

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doa_threshold = 10;

full_test_sim_az(doa_true_angles_az, d_row, T, lambda, Nr, N, SNR, rand_seed, doa_threshold);
% full_test_sim_elev(doa_true_angles_elev, d_col, T, lambda, Nr, N, SNR, rand_seed, doa_threshold);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = full_test_sim_elev(doa_true_angles_elev, d_col, T, lambda, Nr, N, SNR, rand_seed, doa_threshold)
    % 1. Generate simulated measurements
    test_signal_sim_elev = simulate_data(doa_true_angles_elev, T, lambda, d_col, Nr, SNR, N, rand_seed);
    
    % 2. Filter samples in each client & frame which are below given DoA threshold
    % Create array of cells of a size (M x N), where M - # of clients, N - # of frames
    test_signal_sim_elev_filtered = cell(size(test_signal_sim_elev, 3), size(test_signal_sim_elev, 4));
    for client_idx = 1:size(test_signal_sim_elev_filtered, 1)
        fprintf(strcat('\n\nClient: ', int2str(client_idx)));
        for frame_idx = 1:size(test_signal_sim_elev_filtered, 2)
            fprintf(strcat('\nFrame: ', int2str(frame_idx)));
    
            % Determine which samples within this frame are "within threshold"
            client_true_doa = doa_true_angles_elev(client_idx);
            client_samples = squeeze(test_signal_sim_elev(:, :, client_idx, frame_idx));
    
            frame_good_samples = find_frame_good_samples_hankel(client_true_doa, client_samples, Nr, doa_threshold, client_idx, false);
    
            % Extract "good" samples from the dataset & save filtered samples into the cell array
            test_signal_sim_elev_filtered(client_idx, frame_idx) = {client_samples(:, frame_good_samples)};
        end
    end
    
    % 3. Run Hanke vs. MUSIC evaluation on filtered dataset
    doa_rmse_hankel_elev = evaluate_hankel(doa_true_angles_elev, test_signal_sim_elev_filtered, Nr, N);
    doa_rmse_music_elev = evalute_music(doa_true_angles_elev, test_signal_sim_elev_filtered, lambda, Nr);
    
    % Plot Hankel performance
    figure
    hold on;
    plot(1:9, doa_rmse_hankel_elev, 'red-o','MarkerFaceColor', 'red', 'DisplayName','Hankel SVD');
    plot(1:9, doa_rmse_music_elev, 'green-o','MarkerFaceColor', 'green', 'DisplayName','MUSIC');
    hold off;
    grid on; 
    xlabel('Client');
    ylabel('DOA estimate RMSE');
    legend();
    title("Hankel vs. MUSIC. Elevation. Simulated Data.");
end

function [] = full_test_sim_az(doa_true_angles_az, d_row, T, lambda, Nr, N, SNR, rand_seed, doa_threshold)
    % 1. Generate simulated measurements
    test_signal_sim_az = simulate_data(doa_true_angles_az, T, lambda, d_row, Nr, SNR, N, rand_seed);
    
    % 2. Filter samples in each client & frame which are below given DoA threshold
    % Create array of cells of a size (M x N), where M - # of clients, N - # of frames
    test_signal_sim_az_filtered = cell(size(test_signal_sim_az, 3), size(test_signal_sim_az, 4));
    for client_idx = 1:size(test_signal_sim_az_filtered, 1)
        fprintf(strcat('\n\nClient: ', int2str(client_idx)));
        for frame_idx = 1:size(test_signal_sim_az_filtered, 2)
            fprintf(strcat('\nFrame: ', int2str(frame_idx)));
    
            % Determine which samples within this frame are "within threshold"
            client_true_doa = doa_true_angles_az(client_idx);
            client_samples = squeeze(test_signal_sim_az(:, :, client_idx, frame_idx));
    
            frame_good_samples = find_frame_good_samples_hankel(client_true_doa, client_samples, Nr, doa_threshold, client_idx, false);
    
            % Extract "good" samples from the dataset & save filtered samples into the cell array
            test_signal_sim_az_filtered(client_idx, frame_idx) = {client_samples(:, frame_good_samples)};
        end
    end
    
    % 3. Run Hanke vs. MUSIC evaluation on filtered dataset
    doa_rmse_hankel_az = evaluate_hankel(doa_true_angles_az, test_signal_sim_az_filtered, Nr, N);
    doa_rmse_music_az = evalute_music(doa_true_angles_az, test_signal_sim_az_filtered, lambda, Nr);
    
    % Plot Hankel performance
    figure
    hold on;
    plot(1:9, doa_rmse_hankel_az, 'red-o','MarkerFaceColor', 'red', 'DisplayName','Hankel SVD');
    plot(1:9, doa_rmse_music_az, 'green-o','MarkerFaceColor', 'green', 'DisplayName','MUSIC');
    hold off;
    grid on; 
    xlabel('Client');
    ylabel('DOA estimate RMSE');
    legend();
    title("Hankel vs. MUSIC. Azimuth. Simulated Data.");
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

function [doa_rmse] = evalute_music(doa_true_angles, test_signal_sim, lambda, Nr)
    doa_rmse = zeros(size(test_signal_sim, 1), 1);
    for doa_idx = 1:length(doa_true_angles)

        true_doa = doa_true_angles(doa_idx);

        frame_errors = zeros(size(test_signal_sim, 2), 1);
        for frame_idx = 1:size(test_signal_sim, 2)
            X = cell2mat(test_signal_sim(doa_idx, frame_idx));

            doa_estimated = music.run_music(X, 1, lambda, Nr);
            frame_errors(frame_idx, 1) = (doa_estimated - true_doa)^2;
        end
        doa_rmse(doa_idx, 1) = sqrt(mean(frame_errors));
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

            % data_matrix = test_signal_sim(:, :, doa_idx, frame_idx);
            
            sample_rmse = zeros(size(data_matrix, 2), 1);
            for sample_idx = 1:size(data_matrix, 2)
                fprintf('.');
                X = data_matrix(:,sample_idx);
                [svd_final_aoa, ~, ~, ~, ~] = MSSA_SVD_findOptimalWandK_filtering_onetime_max_magnitude(X, Nr, true_doa, 2, 3, 1, 's', 1);
                sample_rmse(sample_idx, 1) = (svd_final_aoa - true_doa)^2;
            end 
            frame_rmse(frame_idx, 1) = sqrt(mean(sample_rmse));
        end
        doa_rmse(doa_idx, 1) = mean(frame_rmse);
    end 
end

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
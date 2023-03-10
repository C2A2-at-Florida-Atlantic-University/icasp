function [final_aoa,final_series,final_min_mse,windowWithMinMSE,optimalKforWindow] = MSSA_SVD_findOptimalWandK_filtering_onetime_theta_min_mse(noisy_series,clean_series,num_antennas,true_aoa,wLower,wUpper,hankMode,orientation)

%from w and k that make theta^ = theta find and return the seq. produced by
%w andk tham min the mse between the filtered seq and the original(clean) 

%input: 
    %noisy_series: data with noise 
    %clean_series: series before the addition of any kind of noise
    %num_antennas: number of antennas 
    %true_aoa: true angle of arrival  
    %wLower: lower bound of window 
    %wUpper: upper bouind of window 
    %wLower and wUpper give the range of different window sizes that we will
    %examine and find the one that produces the optimal denoising and thus the smallest MSE
    %hankMode = 1 average hankelization 
    %hankMode = 0 median hankelization 
    %orientation: what kind of hankelization to use, d or s
    
%mse holds the best mse (lowest) for every window 
%karray holds the k the produced the lowest mse for every window  
%the first element in the mse and karray corresponds to the first window, the second mse and karray corresponds to the second window
%and etc 
window_kArray = [];
window_filtered_series_array = []; %holds the the best filtered series for every window 
window_aoa_array = []; %holds the best aoa for every window
window_min_mse_x_xhat_array = []; %holds the max magnitude for every window 
window_min_error_theta_array = []; 
windowArray = wLower : wUpper; 

%break the complex space series to its constituent real and imaginary part 
real_imag = [real(noisy_series) imag(noisy_series)]; 

%for window size range 
for window = wLower : wUpper
    
    %hankelize each vector time series using the specified window size
    hankelizedData = vectorHankelization(real_imag,window,orientation);
    s = svd(hankelizedData, 'econ');
    
    %these arrays hold the squared error between estimated aoa and true aoa, the filtered series, the corresponding aoa and the error between x and x^ for
    %every rank for a set window 
    rank_filtered_series_array = []; 
    rank_aoa_array = []; 
    rank_mse_x_xhat_array = []; 
    rank_error_theta_array = []; 
    
    %for every singular value
    for k = 1 : length(s)
        %denoised time series given a specific k
        denoised_series = MSSA_SVD_Denoising(real_imag,window,k,2,orientation,hankMode);
        denoised_series = denoised_series.'; 
        denoised_series = complex(denoised_series(:,1),denoised_series(:,2)); 
        
        %calculate the angle of arrival based on the filtered series 
        [estimated_aoa, ~,~] = aoa(denoised_series,num_antennas);
        
        %for every rank of a window store the squared error, the
        %filtered series, the aoa, and the magnitude
        rank_error_theta_array = [rank_error_theta_array ; immse(estimated_aoa,true_aoa)];
        rank_filtered_series_array = [rank_filtered_series_array ; denoised_series.'];
        rank_aoa_array = [rank_aoa_array ; estimated_aoa];
        rank_mse_x_xhat_array = [rank_mse_x_xhat_array ; mean(abs(clean_series-denoised_series).^2)];
        
    end 
    
    %find the indexes of minimum theta error in the rank_error_theta_array (these
    %indexes are the ranks that produced the minimum error between theta^
    %and theta
    min_error_theta_indexes_rank = find(rank_error_theta_array == min(rank_error_theta_array));
    
    %get the mse between x and xhat that correspond to those indexes and
    %find the min
    min_mse_x_xhat_rank = rank_mse_x_xhat_array(min_error_theta_indexes_rank);
    [min_mse_rank,index_min_mse_x_xhat_rank] = min(min_mse_x_xhat_rank);
    
    window_kArray = [window_kArray ; min_error_theta_indexes_rank(index_min_mse_x_xhat_rank)];
    window_filtered_series_array = [window_filtered_series_array ; rank_filtered_series_array(min_error_theta_indexes_rank(index_min_mse_x_xhat_rank),:)];
    window_aoa_array = [window_aoa_array ; rank_aoa_array(min_error_theta_indexes_rank(index_min_mse_x_xhat_rank))];
    window_min_mse_x_xhat_array = [window_min_mse_x_xhat_array ; min_mse_rank];
    window_min_error_theta_array = [window_min_error_theta_array ; min(rank_error_theta_array)];
end 

%find the indexes of the window/s that produces the min theta error 
min_theta_error_indexes_window = find(window_min_error_theta_array == min(window_min_error_theta_array));

%get the mse between x and xhat that correspond to those indexes and
%find the min
min_mse_x_xhat_window = window_min_mse_x_xhat_array(min_theta_error_indexes_window);
[min_mse_window,index_min_mse_x_xhat_window] = min(min_mse_x_xhat_window);

windowWithMinMSE = windowArray(min_theta_error_indexes_window(index_min_mse_x_xhat_window));
optimalKforWindow = window_kArray(min_theta_error_indexes_window(index_min_mse_x_xhat_window)); 
final_aoa = window_aoa_array(min_theta_error_indexes_window(index_min_mse_x_xhat_window));
final_series = window_filtered_series_array(min_theta_error_indexes_window(index_min_mse_x_xhat_window),:).';
final_min_mse = min_mse_window; 
end


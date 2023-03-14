function [final_aoa,final_series,final_max_magnitude,windowWithMinMSE,optimalKforWindow] = MSSA_SVD_findOptimalWandK_filtering_onetime_max_magnitude(noisy_series,num_antennas,true_aoa,wLower,wUpper,hankMode,orientation,step,d, lambda)

%from w and k that make theta^ = theta find and return the seq. produced by
%w and k that max the magnitude. 

%input: 
    %noisy_series: data with noise 
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
window_max_magnitude_array = []; %holds the max magnitude for every window 
window_min_error_array = []; 
windowArray = wLower : wUpper; 

%break the complex space series to its constituent real and imaginary part 
real_imag = [real(noisy_series) imag(noisy_series)]; 

%for window size range 
for window = wLower : wUpper
    
    %hankelize each vector time series using the specified window size
    hankelizedData = vectorHankelization(real_imag,window,orientation);
    s = svd(hankelizedData, 'econ');
    
    %these arrays hold the squared error between estimated aoa and true aoa, the filtered series, the corresponding aoa and the magnitude for
    %every rank for a set window 
    rank_filtered_series_array = []; 
    rank_aoa_array = []; 
    rank_magnitude_array = [];
    rank_error_array = []; 
    
    %for every singular value until MSE is the lowest 
    for k = 1 : length(s)
        %denoised time series given a specific k
        denoised_series = MSSA_SVD_Denoising(real_imag,window,k,2,orientation,hankMode);
        denoised_series = denoised_series.'; 
        denoised_series = complex(denoised_series(:,1),denoised_series(:,2));
        
        %calculate the angle of arrival based on the filtered series 
        [estimated_aoa, magnitude] = aoa(denoised_series,num_antennas,step, d, lambda);
        
        %for every rank of a window store the squared error, the
        %filtered series, the aoa, and the magnitude
        rank_error_array = [rank_error_array ; immse(estimated_aoa,true_aoa)];
        rank_filtered_series_array = [rank_filtered_series_array ; denoised_series.'];
        rank_aoa_array = [rank_aoa_array ; estimated_aoa];
        rank_magnitude_array = [rank_magnitude_array ; magnitude];
        
    end 
    
    %find the indexes of minimum error in the rank_error_array (these
    %indexes are the ranks that produced the minimum error 
    min_error_indexes_rank = find(rank_error_array == min(rank_error_array));
    
    %get the magnitudes that correspond to those indexes and find the max
    min_error_magnitudes_rank = rank_magnitude_array(min_error_indexes_rank);
    [max_magnitude_rank,index_max_magnitude_rank] = max(min_error_magnitudes_rank);
    
    window_kArray = [window_kArray ; min_error_indexes_rank(index_max_magnitude_rank)];
    window_filtered_series_array = [window_filtered_series_array ; rank_filtered_series_array(min_error_indexes_rank(index_max_magnitude_rank),:)];
    window_aoa_array = [window_aoa_array ; rank_aoa_array(min_error_indexes_rank(index_max_magnitude_rank))];
    window_max_magnitude_array = [window_max_magnitude_array ; max_magnitude_rank];
    window_min_error_array = [window_min_error_array ; min(rank_error_array)];
end 

%find the indexes of the window/s that produces the min error 
min_error_indexes_window = find(window_min_error_array == min(window_min_error_array));

%get the magnitudes that correspond to those indexes and find the max 
min_error_magnitudes_window = window_max_magnitude_array(min_error_indexes_window);   
[max_magnitude_window,index_max_magnitude_window] = max(min_error_magnitudes_window);

windowWithMinMSE = windowArray(min_error_indexes_window(index_max_magnitude_window));
optimalKforWindow = window_kArray(min_error_indexes_window(index_max_magnitude_window)); 
final_aoa = window_aoa_array(min_error_indexes_window(index_max_magnitude_window));
final_series = window_filtered_series_array(min_error_indexes_window(index_max_magnitude_window),:).';
final_max_magnitude = max_magnitude_window; 
end


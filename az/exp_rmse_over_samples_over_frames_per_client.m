%clear all ; clc ; close all ; 
step = 1;
clients = 1:8;
ref = 1:8;
num_antennas = 4; 
 
%figure counter 
j = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CLIENT NODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

client_l1_error = []; 
client_svd_error = [];

%for every client 
for c = clients 
    
    %load the workspace 
    work_space_client = 'client_x'+string(c)+'_for_az.mat'; 
    true_aoa_client = double(load(work_space_client).('output_az')); %convert single to double 
    
    client_frame_l1_error = []; 
    client_frame_svd_error = [];
    
    %for every frame  
    for f = 1:19
        
        client_sample_l1_error = []; 
        client_sample_svd_error = [];
        
        %for every snapshot/sample 
        for s = 1:128
            frame_name = 'output_samples_frame_'+string(f); 
            data_tensor = load(work_space_client).(frame_name); %data matrix is a tensor (rows, cols, frame)
            data_matrix = data_tensor(:,:,s); 

            %look for NaNs (NaN matrix is a boolean matrix, if NaN is present will have 1 to that location otherwise 0).
            NaN_matrix = isnan(data_matrix); 
            %find the rows that contain NaNs 
            [row_with_NaN,~] = find(NaN_matrix == 1);

            %removes this row/rows from data_matrix 
            data_matrix(row_with_NaN,:) = [];

            averaged_row = mean(data_matrix,1).';

            [l1_final_aoa,l1_final_series,l1_final_max_magnitude,l1_windowWithMinMSE,l1_optimalKforWindow] = MSSA_L1_findOptimalWandK_filtering_onetime_max_magnitude(averaged_row,num_antennas,true_aoa_client,2,3,1,'s',step);
            [svd_final_aoa,svd_final_series,svd_final_max_magnitude,svd_windowWithMinMSE,svd_optimalKforWindow] = MSSA_SVD_findOptimalWandK_filtering_onetime_max_magnitude(averaged_row,num_antennas,true_aoa_client,2,3,1,'s',step);

            client_sample_l1_error = [client_sample_l1_error ; (l1_final_aoa - true_aoa_client)^2];
            client_sample_svd_error = [client_sample_svd_error ;  (svd_final_aoa - true_aoa_client)^2];
        end 
        
        client_frame_l1_error = [client_frame_l1_error sqrt(mean(client_sample_l1_error))]; 
        client_frame_svd_error = [client_frame_svd_error sqrt(mean(client_sample_svd_error))];
    end
    client_l1_error = [client_l1_error sqrt(mean(client_frame_l1_error))]; 
    client_svd_error = [client_svd_error sqrt(mean(client_frame_svd_error))];
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%REFERENCE NODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ref_l1_error = []; 
ref_svd_error = [];

%for every reference 
for r = ref 
    
    %load the workspace 
    work_space_ref = 'reference_x'+string(r)+'_for_az.mat'; 
    true_aoa_ref = double(load(work_space_ref).('output_az')); %convert single to double 
    
    ref_frame_l1_error = []; 
    ref_frame_svd_error = [];
    
    %for every frame  
    for f = 1:19
        
        ref_sample_l1_error = []; 
        ref_sample_svd_error = [];
        
        %for every snapshot/sample 
        for s = 1:128
            frame_name = 'output_samples_frame_'+string(f); 
            data_tensor = load(work_space_ref).(frame_name); %data matrix is a tensor (rows, cols, frame)
            data_matrix = data_tensor(:,:,s); 

            %look for NaNs (NaN matrix is a boolean matrix, if NaN is present will have 1 to that location otherwise 0).
            NaN_matrix = isnan(data_matrix); 
            %find the rows that contain NaNs 
            [row_with_NaN,~] = find(NaN_matrix == 1);

            %removes this row/rows from data_matrix 
            data_matrix(row_with_NaN,:) = [];

            averaged_row = mean(data_matrix,1).';

            [l1_final_aoa,l1_final_series,l1_final_max_magnitude,l1_windowWithMinMSE,l1_optimalKforWindow] = MSSA_L1_findOptimalWandK_filtering_onetime_max_magnitude(averaged_row,num_antennas,true_aoa_ref,2,3,1,'s',step);
            [svd_final_aoa,svd_final_series,svd_final_max_magnitude,svd_windowWithMinMSE,svd_optimalKforWindow] = MSSA_SVD_findOptimalWandK_filtering_onetime_max_magnitude(averaged_row,num_antennas,true_aoa_ref,2,3,1,'s',step);

            ref_sample_l1_error = [ref_sample_l1_error ; (l1_final_aoa - true_aoa_ref)^2];
            ref_sample_svd_error = [ref_sample_svd_error ;  (svd_final_aoa - true_aoa_ref)^2];
        end 
        
        ref_frame_l1_error = [ref_frame_l1_error sqrt(mean(ref_sample_l1_error))]; 
        ref_frame_svd_error = [ref_frame_svd_error sqrt(mean(ref_sample_svd_error))];
    end
    ref_l1_error = [ref_l1_error sqrt(mean(ref_frame_l1_error))]; 
    ref_svd_error = [ref_svd_error sqrt(mean(ref_frame_svd_error))];
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%l1 client_ref error 
l1_client_ref_error = [client_l1_error ref_l1_error];
%svd client_ref error
svd_client_ref_error = [client_svd_error ref_svd_error];

x_axis = 1:16;
figure(1); 
plot(x_axis,l1_client_ref_error,'b-x','MarkerFaceColor','b');
grid on; 
hold on;
plot(x_axis,svd_client_ref_error,'r-o','MarkerFaceColor','r');
%xlabel('Client');
xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16])
xticklabels({'Client 1','Client 2','Client 3','Client 4','Client 5','Client 6','Client 7','Client 8','Reference 1','Reference 2','Reference 3','Reference 4','Reference 5','Reference 6','Reference 7','Reference 8'});
ylabel('DOA estimate RMSE');
legend('L1','SVD');

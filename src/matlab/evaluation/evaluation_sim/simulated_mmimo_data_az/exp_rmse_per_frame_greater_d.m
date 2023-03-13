clear all; clc; close all; 

PATH_SOURCE = '../../../../data/';

addpath('../../utilities/');

step = 1;
clients = 1:8;
ref = 1;
num_antennas = 4; 
 
%figure counter 
j = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CLIENT NODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%client_l1_error = []; 
client_svd_error = [];

%for every client 
for c = clients 
    fprintf(strcat('\n\nClient', int2str(c)));

    %load the workspace 
    path_client = 'Location_greater_d'+string(c)+'/'; 
    true_aoa_client = double(load(string(path_client)+'Location'+string(c)+'_TrueTheta.mat').('Location'+string(c)+'_TrueTheta')); %convert single to double 
    
    %client_frame_l1_error = []; 
    client_frame_svd_error = [];
    
    %for every frame  
    for f = 1:19
        fprintf(strcat('\nFrame', int2str(f), ':'));
        
        %client_sample_l1_error = []; 
        client_sample_svd_error = [];
        
        %for every snapshot/sample 
        for s = 1:128
            fprintf('.');
            
            frame_name = 'Location'+string(c)+'_Frame'+string(f); 
            data_matrix = load(string(path_client)+string(frame_name)).(frame_name); 

            %[l1_final_aoa,l1_final_series,l1_final_max_magnitude,l1_windowWithMinMSE,l1_optimalKforWindow] = MSSA_L1_findOptimalWandK_filtering_onetime_max_magnitude(data_matrix(:,s),num_antennas,true_aoa_client,2,3,1,'s',step);
            [svd_final_aoa,svd_final_series,svd_final_max_magnitude,svd_windowWithMinMSE,svd_optimalKforWindow] = MSSA_SVD_findOptimalWandK_filtering_onetime_max_magnitude(data_matrix(:,s),num_antennas,true_aoa_client,2,3,1,'s',step);

            %client_sample_l1_error = [client_sample_l1_error ; (l1_final_aoa - true_aoa_client)^2];
            client_sample_svd_error = [client_sample_svd_error ;  (svd_final_aoa - true_aoa_client)^2];
            
        end 
        %client_frame_l1_error = [client_frame_l1_error sqrt(mean(client_sample_l1_error))]; 
        client_frame_svd_error = [client_frame_svd_error sqrt(mean(client_sample_svd_error))];
    end
    %client_l1_error = [client_l1_error sqrt(mean(client_frame_l1_error))]; 
    client_svd_error = [client_svd_error sqrt(mean(client_frame_svd_error))];
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%REFERENCE NODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%client_l1_error = []; 
ref_svd_error = [];

%for every client 
for r = ref 
    fprintf(strcat('\n\nReference', int2str(r)));
    
    %load the workspace 
    path_ref = 'Reference_greater_d/'; 
    true_aoa_ref = double(load(string(path_ref)+'Ref_TrueTheta').('Ref_TrueTheta')); %convert single to double 
    
    %client_frame_l1_error = []; 
    ref_frame_svd_error = [];
    
    %for every frame  
    for f = 1:19
        fprintf(strcat('\nFrame', int2str(f), ':'));
        
        %client_sample_l1_error = []; 
        ref_sample_svd_error = [];
        
        %for every snapshot/sample 
        for s = 1:128
            fprintf('.');
            
            frame_name = 'Ref_Frame'+string(f); 
            data_matrix = load(string(path_ref)+string(frame_name)).(frame_name); 

            %[l1_final_aoa,l1_final_series,l1_final_max_magnitude,l1_windowWithMinMSE,l1_optimalKforWindow] = MSSA_L1_findOptimalWandK_filtering_onetime_max_magnitude(data_matrix(:,s),num_antennas,true_aoa,2,5,1,'s',step);
            [svd_final_aoa,svd_final_series,svd_final_max_magnitude,svd_windowWithMinMSE,svd_optimalKforWindow] = MSSA_SVD_findOptimalWandK_filtering_onetime_max_magnitude(data_matrix(:,s),num_antennas,true_aoa_ref,2,3,1,'s',step);

            %sample_l1_error = [sample_l1_error ; (l1_final_aoa - true_aoa)^2];
            ref_sample_svd_error = [ref_sample_svd_error ;  (svd_final_aoa - true_aoa_ref)^2];
            
        end 
        %frame_l1_error = [frame_l1_error sqrt(mean(sample_l1_error))]; 
        ref_frame_svd_error = [ref_frame_svd_error sqrt(mean(ref_sample_svd_error))];
    end
    %client_l1_error = [client_l1_error sqrt(mean(frame_l1_error))]; 
    ref_svd_error = [ref_svd_error sqrt(mean(ref_frame_svd_error))];
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%l1 client_ref error 
%l1_client_ref_error = [client_l1_error ref_l1_error];
%svd client_ref error
svd_client_ref_error = [client_svd_error ref_svd_error];

%plot
figure(1); 
x_axis = 1:9;
%plot(clients,client_l1_error,'b-x','MarkerFaceColor','b');
%hold on;
plot(x_axis,svd_client_ref_error,'r-o','MarkerFaceColor','r');
grid on; 
xlabel('Client');
ylabel('DOA estimate RMSE');
legend('svd');
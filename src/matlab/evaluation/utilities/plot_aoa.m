function plot_aoa(space_series,num_antennas,step,true_aoa,pred_aoa,magnitude,title_var)
%plot the magnitudes of different angles of arrivals in dB scale

%input: 
    %space_series
    %num_antennas: the number of antennas in the antenna array
    %title: title of plot
    
%output: 
    %magnitudes
    %plot 

    

    %plot 
    magnitudes = [];
    scanning_angle_vector = -89:step:89; 
    for i= scanning_angle_vector
        magnitudes = [magnitudes abs(space_series'*array_response_vector(num_antennas,i))^2];   
    end
    
    %normalized beam pattern of the system in dB
    magnitudes = 10*log10(magnitudes/num_antennas^2);
   
    plot(scanning_angle_vector,magnitudes,'black');
    hold on; 
    plot(pred_aoa,magnitude,'gx');
    xline(true_aoa,'r--');
    title(title_var);
    xlabel('Angle (degrees)');
    ylabel('Normalized Beam Pattern (dB)');
    grid on;
    legend('Spectrum','Estimated AoA: '+string(pred_aoa), 'True AoA: '+string(true_aoa));
end


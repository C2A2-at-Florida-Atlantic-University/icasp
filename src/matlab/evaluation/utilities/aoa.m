function [angle,max_magnitude,magnitudes] = aoa(space_series,num_antennas,step,d, lambda)
%finds the angle of arrival; that is the max(abs(x'*a(theta)))

%input: 
    %space_series
    %num_antennas: the number of antennas in the antenna array
    %step: step for the scanning vector. if 1 every one degree 
%output: 
    %angle: in degrees 
    %max_magnitude
    %magnitudes
    
    magnitudes = [];
    scanning_angle_vector = -89:step:90; 
    for i= scanning_angle_vector
        magnitudes = [magnitudes abs(space_series'*array_response_vector(num_antennas,i,d, lambda))^2];   
    end
    magnitudes = 10*log10(magnitudes/num_antennas^2); %normalized beam pattern of the system in dB
    max_magnitude = max(magnitudes);
    
    %index that max magnitude occurs 
    index = find(magnitudes == max_magnitude);
    
    if length(index)>1 
        index = max(index);
    end
    
    angle = scanning_angle_vector(index);
    
end


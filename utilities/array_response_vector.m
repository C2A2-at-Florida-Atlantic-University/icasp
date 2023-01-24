function a = array_response_vector(num_antennas, theta)
%d=lc/2
%input: 
    %num_antennas
    %theta: angle theta with respect to the broadside 
    
    %generate the array response vector 
    a = 1;
    for index = 1:num_antennas-1
        a = [a;exp(-1j*index*pi*sind(theta))]; 
    end
end



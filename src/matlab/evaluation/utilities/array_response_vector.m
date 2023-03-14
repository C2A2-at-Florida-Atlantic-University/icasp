function a = array_response_vector(num_antennas, theta, d, lambda)
%d=lc/2
%input: 
    %num_antennas
    %theta: angle theta with respect to the broadside 
    
    %generate the array response vector 
    % a = 1;
    % for index = 1:num_antennas-1
    a = exp(-1j*2*pi*(d)*(0:num_antennas-1)'*sind(theta)/lambda);
        % a = [a;exp(-1j*index*pi*sind(theta))]; 
        % a = [a;exp(-1j*2*pi*(d)*(0:Nr-1)'*sind(theta(i))/lambda)]; 
        % SS = ;
    % end
end



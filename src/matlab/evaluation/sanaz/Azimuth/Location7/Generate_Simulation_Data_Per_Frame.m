clc 
clear
close all
%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Location7_TrueTheta = [53];    %Direction of arrival (Degree)
T = 128;   % Samples
K = length(Location7_TrueTheta); %The number of signal source
f = 3.55e9;
lambda = 0.0844486; % Wavelength (meters)
d_x = 0.07935;      % Receiver's antennas spacing (meters)  
Nr  = 4;          %Number of receiver's antennas 
SNR = 1;           %SNR(dB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = zeros(Nr,K);   %Steering Matrix 
for k=1:K 
    A(:,k) = exp(-1j*2*pi*d_x*sind(Location7_TrueTheta(k))*(0:Nr-1)'/lambda); 
end 
%%%%
Vj = diag(sqrt((10.^(SNR/10))/2));
s = Vj* ( randn(K,T) + 1j*randn(K,T) ); %Amplitude vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Location7_Frame1 = A*s+noise;      %Insert AWGN 
save Location7_Frame1
save Location7_TrueTheta

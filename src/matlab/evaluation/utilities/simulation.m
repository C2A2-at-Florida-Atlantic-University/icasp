classdef simulation
    methods (Static)
%         %% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Location2_TrueTheta = [-6.93];    %Direction of arrival (Degree)
%         T = 128;                           % Samples
%         K = length(Location2_TrueTheta); %The number of signal source
%         f = 3.55e9;
%         lambda = 0.0844486; % Wavelength (meters)
%         d_x = 0.07935;      % Receiver's antennas spacing (meters)  0.07935lambda/2;%
%         Nr  = 4;          %Number of receiver's antennas 
%         SNR = 1;           %SNR(dB)
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Simulates MIMO signal given true DoA angle.
        % true_doa: array of DoA values (e.g., [-45])
        % T:        number of samples to generate (e.g., 128)
        % lambda:   wavelength (in meters) (e.g., 0.0844486)
        % d_x:      receiver antenna spacing (in meters) (e.g., 0.07935)
        % Nr:       number of receiver antennas (e.g., 4)
        % SNR:      signal to noise ratio, in dB (e.g., 1)
        function [signal_sim] = simulate_samples(true_doa, T, lambda, d_x, Nr, SNR)
            K = length(true_doa); % determine number of signal sources
            A = zeros(Nr, K); % steering Matrix 
            for k=1:K 
                A(:,k) = exp(-1j*2*pi*d_x*sind(true_doa(k))*(0:Nr-1)'/lambda); 
            end 

            Vj = diag(sqrt((10.^(SNR/10))/2));
            s = Vj * (randn(K,T) + 1j*randn(K,T)); % amplitude vector
            noise = sqrt(1/2) * (randn(Nr,T)+1j*randn(Nr,T));
            
            signal_sim = A * s + noise; % Generate AWGN 
            % save Location2_Frame1
            % save Location2_TrueTheta
        end
    end
end
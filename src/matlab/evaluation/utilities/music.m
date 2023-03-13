classdef music
    methods (Static)
        function [result_doa] = run_music(X, K, lambda, Nr)
            T = size(X, 2);

            Rx = 1 / T * (X * X');

            [eigenVec, ~] = eig(Rx); % find the eigenvalues and eigenvectors of Rx 
            Vn = eigenVec(:,1:Nr-K);          
            theta = -90:1:90;              
            for i=1:length(theta) 
                SS = zeros(Nr,1); 
                SS = exp(-1j*2*pi*(lambda/2)*(0:Nr-1)'*sind(theta(i))/lambda);
                Pmusic(i) = 1/(SS'*(Vn*Vn')*SS); 
            end
            Pmusic = real(10*log10(Pmusic)); % spatial spectrum function
            [~, locs] = findpeaks(Pmusic, theta, 'SortStr', 'descend', 'Annotate', 'extents');

            if isempty(locs)
                result_doa = NaN;
            else
                result_doa = locs(1);
            end
        end
    end
end
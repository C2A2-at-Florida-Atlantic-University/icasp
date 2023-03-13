clc 
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load 'Location8_TrueTheta.mat', load 'Location8_Frame1.mat', load 'Location8_Frame2.mat', load 'Location8_Frame3.mat',  load 'Location8_Frame4.mat', load 'Location8_Frame5.mat',load 'Location8_Frame6.mat', load 'Location8_Frame7.mat', load 'Location8_Frame8.mat', load 'Location8_Frame9.mat', load 'Location8_Frame10.mat', load 'Location8_Frame11.mat', load 'Location8_Frame12.mat', load 'Location8_Frame13.mat', load 'Location8_Frame14.mat', load 'Location8_Frame15.mat', load 'Location8_Frame16.mat', load 'Location8_Frame17.mat', load 'Location8_Frame18.mat', load 'Location8_Frame19.mat'
DataFrames = [{Location8_Frame1} {Location8_Frame2} {Location8_Frame3} {Location8_Frame4} {Location8_Frame5} {Location8_Frame6} {Location8_Frame7} {Location8_Frame8} {Location8_Frame9} {Location8_Frame10} {Location8_Frame11} {Location8_Frame12} {Location8_Frame13} {Location8_Frame14} {Location8_Frame15} {Location8_Frame16} {Location8_Frame17} {Location8_Frame18} {Location8_Frame19}];
%load the workspace
  for f = 1:19
    true_aoa = Location8_TrueTheta;
    X = DataFrames{f};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MUSIC 
Rx = cov(X');                     %Data covarivance matrix 
[eigenVec,eigenVal] = eig(Rx);    %Find the eigenvalues and eigenvectors of Rx 
Vn = eigenVec(:,1:Nr-K);          
theta = -90:1:90;              
for i=1:length(theta) 
    SS = zeros(Nr,1); 
    SS = exp(-1j*2*pi*(lambda/2)*(0:Nr-1)'*sind(theta(i))/lambda);
    Pmusic(i) = 1/(SS'*(Vn*Vn')*SS); 
end
Pmusic = real(10*log10(Pmusic)); %Spatial Spectrum function
[pks,locs] = findpeaks(Pmusic,theta,'SortStr','descend','Annotate','extents');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
estimated_Angel = locs(1);
error_per_frame(f) = (norm(estimated_Angel-true_aoa))^2
  end
RMSE8 = sqrt(mean(error_per_frame))
save RMSE8



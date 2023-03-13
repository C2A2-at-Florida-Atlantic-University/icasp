clc 
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load 'Location5_TrueTheta.mat', load 'Location5_Frame1.mat', load 'Location5_Frame2.mat', load 'Location5_Frame3.mat',  load 'Location5_Frame4.mat', load 'Location5_Frame5.mat',load 'Location5_Frame6.mat', load 'Location5_Frame7.mat', load 'Location5_Frame8.mat', load 'Location5_Frame9.mat', load 'Location5_Frame10.mat', load 'Location5_Frame11.mat', load 'Location5_Frame12.mat', load 'Location5_Frame13.mat', load 'Location5_Frame14.mat', load 'Location5_Frame15.mat', load 'Location5_Frame16.mat', load 'Location5_Frame17.mat', load 'Location5_Frame18.mat', load 'Location5_Frame19.mat'
DataFrames = [{Location5_Frame1} {Location5_Frame2} {Location5_Frame3} {Location5_Frame4} {Location5_Frame5} {Location5_Frame6} {Location5_Frame7} {Location5_Frame8} {Location5_Frame9} {Location5_Frame10} {Location5_Frame11} {Location5_Frame12} {Location5_Frame13} {Location5_Frame14} {Location5_Frame15} {Location5_Frame16} {Location5_Frame17} {Location5_Frame18} {Location5_Frame19}];
%load the workspace
  for f = 1:19
    true_aoa = Location5_TrueTheta;
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
RMSE5 = sqrt(mean(error_per_frame))
save RMSE5
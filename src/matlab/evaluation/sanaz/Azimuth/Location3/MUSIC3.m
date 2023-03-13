clc 
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load 'Location3_TrueTheta.mat', load 'Location3_Frame1.mat', load 'Location3_Frame2.mat', load 'Location3_Frame3.mat',  load 'Location3_Frame4.mat', load 'Location3_Frame5.mat',load 'Location3_Frame6.mat', load 'Location3_Frame7.mat', load 'Location3_Frame8.mat', load 'Location3_Frame9.mat', load 'Location3_Frame10.mat', load 'Location3_Frame11.mat', load 'Location3_Frame12.mat', load 'Location3_Frame13.mat', load 'Location3_Frame14.mat', load 'Location3_Frame15.mat', load 'Location3_Frame16.mat', load 'Location3_Frame17.mat', load 'Location3_Frame18.mat', load 'Location3_Frame19.mat'
DataFrames = [{Location3_Frame1} {Location3_Frame2} {Location3_Frame3} {Location3_Frame4} {Location3_Frame5} {Location3_Frame6} {Location3_Frame7} {Location3_Frame8} {Location3_Frame9} {Location3_Frame10} {Location3_Frame11} {Location3_Frame12} {Location3_Frame13} {Location3_Frame14} {Location3_Frame15} {Location3_Frame16} {Location3_Frame17} {Location3_Frame18} {Location3_Frame19}];
%load the workspace
  for f = [1:4 6:19]
    true_aoa = Location3_TrueTheta;
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
RMSE3 = sqrt(mean(error_per_frame))
save RMSE3
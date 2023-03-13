clc 
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load 'Location2_TrueTheta.mat', load 'Location2_Frame1.mat', load 'Location2_Frame2.mat', load 'Location2_Frame3.mat',  load 'Location2_Frame4.mat', load 'Location2_Frame5.mat',load 'Location2_Frame6.mat', load 'Location2_Frame7.mat', load 'Location2_Frame8.mat', load 'Location2_Frame9.mat', load 'Location2_Frame10.mat', load 'Location2_Frame11.mat', load 'Location2_Frame12.mat', load 'Location2_Frame13.mat', load 'Location2_Frame14.mat', load 'Location2_Frame15.mat', load 'Location2_Frame16.mat', load 'Location2_Frame17.mat', load 'Location2_Frame18.mat', load 'Location2_Frame19.mat'
DataFrames = [{Location2_Frame1} {Location2_Frame2} {Location2_Frame3} {Location2_Frame4} {Location2_Frame5} {Location2_Frame6} {Location2_Frame7} {Location2_Frame8} {Location2_Frame9} {Location2_Frame10} {Location2_Frame11} {Location2_Frame12} {Location2_Frame13} {Location2_Frame14} {Location2_Frame15} {Location2_Frame16} {Location2_Frame17} {Location2_Frame18} {Location2_Frame19}];
%load the workspace
  for f = 1:19
    true_aoa = Location2_TrueTheta;
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
RMSE2 = sqrt(mean(error_per_frame))
save RMSE2



clc 
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load 'Ref_TrueTheta.mat', load 'Ref_Frame1.mat', load 'Ref_Frame2.mat', load 'Ref_Frame3.mat',  load 'Ref_Frame4.mat', load 'Ref_Frame5.mat',load 'Ref_Frame6.mat', load 'Ref_Frame7.mat', load 'Ref_Frame8.mat', load 'Ref_Frame9.mat', load 'Ref_Frame10.mat', load 'Ref_Frame11.mat', load 'Ref_Frame12.mat', load 'Ref_Frame13.mat', load 'Ref_Frame14.mat', load 'Ref_Frame15.mat', load 'Ref_Frame16.mat', load 'Ref_Frame17.mat', load 'Ref_Frame18.mat', load 'Ref_Frame19.mat'
DataFrames = [{Ref_Frame1} {Ref_Frame2} {Ref_Frame3} {Ref_Frame4} {Ref_Frame5} {Ref_Frame6} {Ref_Frame7} {Ref_Frame8} {Ref_Frame9} {Ref_Frame10} {Ref_Frame11} {Ref_Frame12} {Ref_Frame13} {Ref_Frame14} {Ref_Frame15} {Ref_Frame16} {Ref_Frame17} {Ref_Frame18} {Ref_Frame19}];
%load the workspace
  for f = 1:19
    true_aoa = Ref_TrueTheta;
    X = DataFrames{f};
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
estimated_Angel = locs(1);
error_per_frame(f) = (norm(estimated_Angel-true_aoa))^2
  end
RMSERef = sqrt(mean(error_per_frame))
save RMSERef



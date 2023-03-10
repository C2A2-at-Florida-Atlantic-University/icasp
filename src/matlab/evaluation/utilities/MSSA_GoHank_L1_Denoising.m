function finalFeatures = MSSA_GoHank_L1_Denoising(data,window,k,numFeatures,orientation,mode)
%input: 
    %data: a matrix where each row is an observation and each col is a
    %feature 
    %window:window 
    %k: optimal k 
    %numFeatures: number of cols if data matrix has been preprocessed 
    %orientation: what kind of hankelization is being used, d or s 
    
%output: 
    %a matrix where each row is the sequence for the corresponding feature
    
    %sample support 
    ss = size(data,1);
    
    %hankelize each vector time series using the specified window size
    hankelizedData = vectorHankelization(data,window,orientation);
    
    %apply the L1 Uniform Feature Preservation
    [Q,~,~,Z] = L1PCA_uniform_feature_preservation(hankelizedData,k);
    
    %hankelization 
    if (mode == 1) %avg hankelization 
        hankelizedLowerRankApproximation = vectorAvgHankelize(Q*Z,window,numFeatures,ss,orientation);
    else 
        hankelizedLowerRankApproximation = vectorMedianHankelize(Q*Z,window,numFeatures,ss,orientation);
    end 
    
    %take out the sequence for each feature
    %output: a matrix were each row is the sequence for the
    %corresponding feature
    finalFeatures = takeFeaturesOut(hankelizedLowerRankApproximation,ss,window,numFeatures,orientation); 
end


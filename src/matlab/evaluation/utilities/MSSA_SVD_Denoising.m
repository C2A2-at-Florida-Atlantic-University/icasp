function finalFeatures = MSSA_SVD_Denoising(data,window,k,numFeatures,orientation,mode)
%input: 
    %data: a matrix where each row is an observation and each col is a
    %feature 
    %window:window 
    %k: optimal k 
    %numFeatures: number of cols if data matrix has been preprocessed 
    %orientation: what kind of hankelization is being used, d or s 
    %mode: 1 for avg hankelization, 0 for median hankelization 
    
%output: 
    %a matrix where each row is the sequence for the corresponding feature
    
    %sample support 
    ss = size(data,1);
    
    %hankelize each vector time series using the specified window size
    hankelizedData = vectorHankelization(data,window,orientation);
    
    %compute the SVD
    %NOTE: svd function will return the USV NOT the USV'
    [U,S,V] = svd(hankelizedData,'econ');

    approximationOfHankelizedY = U(:,1:k)*S(1:k,1:k)*V(:,1:k)';
    
    %hankelization 
    if (mode == 1) 
        hankelizedLowerRankApproximation = vectorAvgHankelize(approximationOfHankelizedY,window,numFeatures,ss,orientation);
    else 
        hankelizedLowerRankApproximation = vectorMedianHankelize(approximationOfHankelizedY,window,numFeatures,ss,orientation);
    end 
    
    %take out the sequence for each feature
    %output: a matrix were each row is the sequence for the
    %corresponding feature
    finalFeatures = takeFeaturesOut(hankelizedLowerRankApproximation,ss,window,numFeatures,orientation); 
end


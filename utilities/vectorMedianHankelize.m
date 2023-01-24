function medianMatrix = vectorMedianHankelize(productMatrix,window,numFeatures,ss,orientation)
%this function will get as input a matrix that has the same size with the
%data representation matrix that was originally created from the dataset. This matrix in
%particular will be the USV product. NOTE: the USV product has the same
%dimensions with the original data representation matrix. 


% GOAL: Taking the USV product matrix and knowing the window size that
% was used to create the data representation matrix as well as the
% numFeatures MEDIAN HANKELIZE each feature individually. 
 
 %avgMatrix = []; 
 medianMatrix = zeros(size(productMatrix,1),size(productMatrix,2)); 
    
    if orientation == 'd'
        %take the features out of the data representation matrix
        %jump will account for the rows we need to skip to get each feature out
        %of the matrix. For the first feature we do not have to jump rows for
        %the second we need to jump window rows for the third 2 windows 
        jump = 0; 
        for feature = 1:numFeatures

            featureMatrix = productMatrix(1+jump:window*feature,:);
            medianHankelized = medianHankelize(featureMatrix);
            medianMatrix(1+jump:window*feature,:) = medianHankelized; 
            %avgMatrix = [avgMatrix ; avgHankelized];

            %update the jump
            jump = jump + window ; 
        end 
    end
    
    if orientation == 's'
        jump = 0;
        for feature = 1:numFeatures

            featureMatrix = productMatrix(:,1+jump:feature*(ss - window +1));
            medianHankelized = medianHankelize(featureMatrix);
            %avgMatrix = [avgMatrix avgHankelized];
            medianMatrix(:,1+jump:feature*(ss - window +1)) = medianHankelized; 
            
            %update the jump
            %in the created matrix each feature spans N-WINDOW +1 cols
            jump = jump + (ss - window +1) ;
        end
    end 
end 


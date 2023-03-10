function  featuresMatrix = takeFeaturesOut(matrix,ss,window,numFeatures,orientation)
    %input: matrix after USVproduct has been calculated
    %ss: sample support available for each coordinate level scalar time
    %series 
    %output: take the sequence for each feature out of the hankel matrix. The output will be a matrix where 
    %each row vector in the matrix will represent the sequence for the corresponding feature. 
    %NOTE: the sequence for each feature will come out in the same order that
    %the vectorHankelization was created
    
    featuresMatrix = []; 
    
    if orientation == 'd'
        %take the features out of the hankel matrix
        %jump will account for the rows we need to skip to get each feature out
        %of the matrix. For the first feature we do not have to jump rows for
        %the second we need to jump window rows for the third 2 windows 
        jump = 0; 
        for feature = 1:numFeatures
            %read the first col up to the window size 
            column = matrix(1+jump:window*feature, 1);

            %read the last row. Do NOT read the first element of the last row because
            %it is already in the output. It is in the first column 
            row = matrix(window+jump,2:end);

            %construct the feature 
            individualFeature = [column.' row];

            %update the jump
            jump = jump + window ; 

            featuresMatrix = [featuresMatrix ; individualFeature];
        end 
    end
    
    if orientation == 's'
        %take the features out of the hankel matrix
        %jump will account for the number of cols we need to skip to get each feature out
        %of the matrix. For the first feature we do not have to jump cols for
        %the second we need to jump N-WINDOW+1 cols. 
        jump = 0;
        for feature = 1:numFeatures
            
            %read the cols
            column = matrix(1:size(matrix,1),1+jump);

            %read the last row. Do NOT read the first element of the last row because
            %it is already in the output. It is in the first column 
            row = matrix(end,2+jump:feature*(ss - window +1));

            %construct the feature 
            individualFeature = [column.' row];

            %update the jump
            %in the created matrix each feature spans N-WINDOW +1 cols
            jump = jump + (ss - window +1) ; 

            featuresMatrix = [featuresMatrix ; individualFeature];
        end 
    end 
end


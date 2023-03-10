function dataRepresentation = vectorHankelization(data,window,orientation)
    %input: data is a matrix with N observations and p features for each
    %observation
    %orientation: d=down stack on top of each other. s = side put side by
    %side 
    %output: a matrix of size (p*D)x(N-D+1). Hankelize each col (feature) with a specific window D
    %and stack all the hankelized features in a matrix if orientation is
    %set to 'd' else if orientation is set to 's' put them next to each
    %other 

    %this will hold all the individual hankelized features 
    dataRepresentation = [];
    
    if orientation == 'd'
        %for every col(feature) in the dataset
        for col = 1:size(data,2)
            %for every feature (col) create a hankel matrix with a specific
            %window
            hankelizedFeature = hankelMatrix(data(:,col),window);
            dataRepresentation = [dataRepresentation ; hankelizedFeature]; 
        end   
    end
    
    if orientation == 's'
        %for every col(feature) in the dataset
        for col = 1:size(data,2)
            %for every feature (col) create a hankel matrix with a specific
            %window
            hankelizedFeature = hankelMatrix(data(:,col),window);
            dataRepresentation = [dataRepresentation hankelizedFeature]; 
        end   
    end
end


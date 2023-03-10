function hankelizedMatrix = avgHankelize(matrix)
%hankelize the given vector by filling the anti-diagonals with the avg
%value of the corresponding anti-diagonal of the given matrix
%we assume that the window does not change and can not be specified 

%flip the matrix
flipedMatrix = flip(matrix);

%split the size of the matrix to rows x cols 
%Diagonal number, specified as an integer. k=0 represents the main diagonal, k>0 is above the main diagonal, and k<0 is below the main diagonal.
%For an m-by-n matrix, k is in the range (−m+1)≤k≤(n−1)  .

numRows = size(matrix,1); 
numCols = size(matrix,2);

%row and column arrays that will hold the corresponding avgs to compute the hankel matrix hankel(r,c)
c = zeros(size(matrix,1),0);
r = zeros(size(matrix,2),0);

%since the matrix is flipped the diag function will return the
%anti-diagonal
%compute c 
for col = (-numRows+1):0
    antidiagonal = diag(flipedMatrix, col);
    c(abs(col)+1)= mean(antidiagonal);
end

%flip c to have the right order of elements in the array 
c = flip(c);

%compute the r 
for row = 0:(numCols-1)
    antidiagonal = diag(flipedMatrix, row); 
    r(row+1)= mean(antidiagonal);
end

%pass it to the hankel function to get the hankelized matrix 
hankelizedMatrix = hankel(c,r);

end


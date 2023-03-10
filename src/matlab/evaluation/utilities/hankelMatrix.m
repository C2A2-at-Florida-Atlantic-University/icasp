function hankelMatrix = hankelMatrix(input,windowSize)
%takes an input array or row vector and a window size and creates a hankel matrix
r = input(windowSize:end);
c = input(1:windowSize);
hankelMatrix = hankel(c,r);
end


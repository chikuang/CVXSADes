function We = initializeExact2(W, nn)
% initializeExact(W,N): given list of weights W for approximate design, and
% total number of design points N, return modified weight list We such that
% for each w[i] in We, N*w gives the exact integer number of points for
% the i-th support point in the exact design
 
We = floor(nn*W); % take floor of approx design*N, sum(We) <= N initially
idx = We==0;
out=sum(idx(:));
for i = 1:length(We)
  if idx(i) == 1
    [M,I] = max( floor(We));
    We(I) = We(I) - 1;
    We(i) = We(i) + 1;
  end
end

i = 1;
while(sum(We) < nn)  
    We(i) = We(i) + 1; % lazily increase number of         
    i = i + 1;         % support points until we have N of them
end

We = We/nn;

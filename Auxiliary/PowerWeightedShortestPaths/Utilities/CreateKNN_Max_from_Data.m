function A = CreateKNN_Max_from_Data(Data,k,r)
%        A = CreateKNN_Max_from_Data(X,K,r)
% Create the adjacency matrix of a (weighted) K-Nearest Neighbours graph
% from a DATA MATRIX X. Gaussian kernel is used.
% Local scaling, based on `Self-Tuning for Spectral Clustering' by
% Zelnick-Manor and Perona, is also used.
% Symmetrize by taking a max.
% 
% INPUT
% =========================================
% Data ................. n-by-d data matrix. Data points stored as rows
% K .................... Number of nearest neighbours to keep.
% r .................... Local clustering parameters will be set using the
% r-th nearest neighbour
%
%OUTPUT
% =========================================
% A ................... Weighted adjacency matrix.
%
% Daniel Mckenzie, using some code from Ke Yin's graph_affinity_matrix.m
% 18th December 2018

[n,d] = size(Data);
[IDX,D] = knnsearch(Data,Data,'K',k,'NSMethod','kdtree');
Scales = D(:,r);  % find local scaling parameters
Dists = reshape(D',[],1);
I = reshape(ones(k, 1) * (1:n), [], 1);
J = reshape(IDX',[],1);
Sigmas = Scales(I).*Scales(J);
S = exp(-Dists.^2./Sigmas);

Atemp = sparse(I,J,S,n,n);
A = max(Atemp',Atemp);

end

    
    
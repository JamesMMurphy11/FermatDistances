function A = NaiveShortestPathAdjMat(X,k,r,p)
%        A = NaiveShortestPathAdjMat(X,k,r,p)
% Inefficient implementation of finding the k path-nearest neighbors.
% Used to verify that ShortestPathAdjMat is working.

[n,d] = size(X);
[IDX,D] = knnsearch(X,X,'K',k,'NSMethod','kdtree');
D = D.^p;  % power weight the distances

% ======== Make a sparse matrix to represent (directed) graph ======== %

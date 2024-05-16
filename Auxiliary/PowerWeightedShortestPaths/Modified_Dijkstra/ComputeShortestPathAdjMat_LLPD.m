function A = ComputeShortestPathAdjMat_LLPD(X,k,r,sym)
%        A = ComputeShortestPathAdjMat_LLPD(X,k,r,sym)
% This function creates a symmetrized, k nearest neighbor matrix using
% Zelnik-Manor local scaling, for an input parameter k. Note that the 
% nearest neighbors are calculated using the longest leg path distance.
% Provides two possible symmetrizations, AMAx and AMult
%
% INPUT
% ===============================================
% X .................. n-by-d matrix. Data points in R^d stored as rows.
% k .................. number of nearest neighbors to keep.
% r .................. local scaling parameter
% sym ................ Choose the symmetrization, either 'max' or 'mult'
%
% OUTPUT
% ==============================================
% A .................. A LLPD nearest neighbors adjacency matrix,
% symmetrized according to sym
%
% Daniel Mckenzie 
% 24 January 2019

[n,d] = size(X);
[IDX,D] = knnsearch(X,X,'K',k,'NSMethod','kdtree');

% ============ Preallocate to make the sparse matrix A ========== %
I = zeros(n*k,1);
J = zeros(n*k,1);
S_temp = zeros(n*k,1);
Scales = zeros(n,1); %for ZMP local scaling

% ============= Loop through all vertices, find their k nearest neighbors
% in path distance ================= %
for i = 1:n
    [K_nearest_neighbors, Distances] = Dijkstra_with_early_stopping_LLPD(IDX,D,i);
    I((i-1)*k + 1:i*k) = i;
    J((i-1)*k + 1:i*k) = K_nearest_neighbors;
    S_temp((i-1)*k + 1:i*k) = Distances;
    Scales(i) = Distances(r);
end

Sigmas = Scales(I).*Scales(J);
S = exp(-S_temp.^2./Sigmas);
Atemp = sparse(I,J,S,n,n);

if sym == 'max'
    A = max(Atemp,Atemp');
elseif sym == 'mult'
    A = Atemp'*Atemp;
else
    error('Error: Choose a symmetrization')
end

end

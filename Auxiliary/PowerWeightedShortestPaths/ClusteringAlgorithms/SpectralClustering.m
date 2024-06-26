function [C, L, U, Eigvals, W] = SpectralClustering(W, k, Type, SpectralOpts)
%SPECTRALCLUSTERING Executes spectral clustering algorithm
%   Executes the spectral clustering algorithm defined by
%   Type on the adjacency matrix W and returns the k cluster
%   indicator vectors as columns in C.
%   If L and U are also called, the (normalized) Laplacian and
%   eigenvectors will also be returned.
%
%   'W' - Adjacency matrix, needs to be square
%   'k' - Number of clusters to look for
%   'Type' - Defines the type of spectral clustering algorithm
%            that should be used. Choices are:
%      1 - Unnormalized
%      2 - Normalized according to Shi and Malik (2000)
%      3 - Normalized according to Jordan and Weiss (2002)
%      4 - Normalized according to Jordan and Weiss AND regularized 
%
%   References:
%   - Ulrike von Luxburg, "A Tutorial on Spectral Clustering", 
%     Statistics and Computing 17 (4), 2007
%
%   Author: Ingo Buerk
%   Year  : 2011/2012
%   Bachelor Thesis

% calculate degree matrix
degs = sum(W, 2);
D    = sparse(1:size(W, 1), 1:size(W, 2), degs);

if strcmp(SpectralOpts.DMNormalization,'yes')
    N = spdiags(1./(degs.^SpectralOpts.DMalpha), 0, size(D, 1), size(D, 2));
    W = N*W*N;
    % Recalculate degrees using new kernel:
    degs = sum(W, 2);
    D    = sparse(1:size(W, 1), 1:size(W, 2), degs);
end

% compute unnormalized Laplacian
L = D - W;

% compute normalized Laplacian if needed
switch Type
    case 2
        % avoid dividing by zero
        degs(degs == 0) = eps;
        % calculate inverse of D
        D = spdiags(1./degs, 0, size(D, 1), size(D, 2));
        
        % calculate normalized Laplacian
        L = D * L;
    case 3
        % avoid dividing by zero
        degs(degs == 0) = eps;
        % calculate D^(-1/2)
        D = spdiags(1./(degs.^0.5), 0, size(D, 1), size(D, 2));
        
        % calculate normalized Laplacian
        L = D * L * D;
    case 4
        % avoid dividing by zero
        degs(degs == 0) = eps;
        % calculate D = (Deg +\tau*I)^(-1/2)
        tau = mean(degs)/10;
        D = spdiags(1./((degs+tau).^0.5), 0, size(D, 1), size(D, 2));
        
        % calculate normalized Laplacian
        L = eye(size(W,1))- D * W * D; 
end

% compute the eigenvectors corresponding to the k smallest
% eigenvalues
diff   = eps;
[U, Eigvals] = eigs(L, k, diff);

% in case of the Jordan-Weiss algorithm, we need to normalize
% the eigenvectors row-wise
if Type == 3
    U = bsxfun(@rdivide, U, sqrt(sum(U.^2, 2)));
end

% now use the k-means algorithm to cluster U row-wise
% C will be a n-by-1 matrix containing the cluster number for
% each data point
%C = kmeans(U, k, 'start', 'cluster', ...
%                 'EmptyAction', 'singleton','Replicates',10,'MaxIter',1000);
C = kmeans(U, k,'Replicates',30);
             
% now convert C to a n-by-k matrix containing the k indicator
% vectors as columns
size_D_1 = size(D, 1);
C = sparse(1:size(D, 1), C, 1);

end
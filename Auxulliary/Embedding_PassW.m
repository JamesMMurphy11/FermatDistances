function [Vecs,Vals,Labels,K_Estimate_Gap,W,L,time,sigma] = Embedding_PassW(DistMatrix,epsilon,NumEigs,K,KernelName)

time=cputime;

sigma = prctile(DistMatrix(DistMatrix>0),100*epsilon);

if strcmp(KernelName,'Gaussian')
    W = exp(-DistMatrix.^2./(sigma^2)).*(DistMatrix>0);
elseif strcmp(KernelName, 'Cutoff')
    W = (DistMatrix<sigma).*(DistMatrix>0);
end

W = max(W,W'); %symmetrize

rootD = diag(sum(W, 2).^(-.5));
L = eye(size(W))-rootD*W*rootD;
% Compute eigendecomposition, normalize, cluster, and get accuracy
%[U, Vals] = eigs(L, NumEigs, eps);
[U, Vals] = eigs(L, NumEigs,'smallestreal');
Vecs = bsxfun(@rdivide, U, sqrt(sum(U.^2, 2)));
Vals=diag(Vals);
Vals=sort(Vals);
[~,K_Estimate_Gap]=max(diff(Vals(2:end)));
K_Estimate_Gap=K_Estimate_Gap+1;

Labels=kmeans(real(Vecs(:,1:K)),K,'Replicates',100);

time = cputime-time;
end


function [FD, ED_Rescaled, ED] = GetSpectra(X,p,NumEig)

%% Compute Laplacians with the cutoff kernel

m=p+1;
q=2*p;
r=0;
n=size(X,1);

tic;
[EucNN_Indices, ED_kNN] = knnsearch(X,X,'k',n);
EucNN_Distances_p = ED_kNN.^p;
Time_Euc_NN=toc;

%% Set a lower bound on the connectivity radius to ensure the Euclidean graph is well-connected

h_ED=max(ED_kNN(:,floor(log(n))));

%% Compute Fermat distances

FD_idx = zeros(size(EucNN_Indices));
FD_kNN = inf*ones(n,n);

NeedsNeighbors=ones(size(X,1),1);
tic;

kNN_Dijkstra=30*floor(log(n));  % Start with this
Count=0;
h_FD=inf;

while any(NeedsNeighbors)
    NumberSearchesThisRound=sum(NeedsNeighbors);
    Count=Count+1;
    for i=1:n
        if NeedsNeighbors(i)
            [temp_idx, temp] = Dijkstra_with_early_stopping(EucNN_Indices(:,1:kNN_Dijkstra), EucNN_Distances_p(:,1:kNN_Dijkstra), i);
            FD_idx(i,1:kNN_Dijkstra) = temp_idx';
            FD_kNN(i,1:kNN_Dijkstra) = temp';
            display(['Round ', num2str(Count), ', Computing Fermat Distance for node ', num2str(i), ' total number of searches needed is ', num2str(NumberSearchesThisRound)]);
        end
        if Count>1
            if FD_kNN(i,kNN_Dijkstra)>h_FD
                NeedsNeighbors(i)=0;
            end
        end
    end

    if Count==1  % In the first pass, get the h value
        h_FD=max(FD_kNN(:,floor(log(n))));
        for i=1:n
            if FD_kNN(i,kNN_Dijkstra)>h_FD
                NeedsNeighbors(i)=0;
            end
        end
    end

    kNN_Dijkstra=min(2*kNN_Dijkstra,n);
end

Time_FD_NN=toc;

Base=[];
for j=1:n
    Base=vertcat(Base,j*ones(sum(FD_kNN(j,:)<Inf),1));
end

FD_KNN_T=FD_kNN';
FD_IDX_T=FD_idx';
FD_IDX_T=FD_IDX_T(FD_KNN_T<Inf);
FD_KNN_T=FD_KNN_T(FD_KNN_T<Inf); % Remove space holders


%% Compute FD Laplacian

tic;

FD_WeightMatrix = sparse(Base(:),FD_IDX_T(:),n^(-1)*h_FD^(-2).*(FD_KNN_T(:)<h_FD),n,n);
FD_WeightMatrix = max(FD_WeightMatrix, FD_WeightMatrix');
FD_WeightMatrix(logical(eye(size(FD_WeightMatrix))))=0;
FD_DegreeMatrix = diag(sum(FD_WeightMatrix,2));

if min(diag(FD_DegreeMatrix))==0
    pause
end

L_FD = eye(size(FD_WeightMatrix))-FD_DegreeMatrix^(-1)*FD_WeightMatrix; % 2-3-23: RW normalization

Time_Build_FD_L=toc;

%% Compute Euclidean Weight Matrix and Laplacian

tic;

ED_WeightMatrix = sparse((n)^(-1)*h_ED^(-2).*(squareform(pdist(X))<h_ED));
ED_WeightMatrix = max(ED_WeightMatrix,ED_WeightMatrix');
ED_WeightMatrix(logical(eye(size(ED_WeightMatrix))))=0;
ED_Degrees = sum(ED_WeightMatrix,2);

if min(ED_Degrees)==0
    pause;
end

ED_DegreeMatrix=diag(ED_Degrees);
L_ED=eye(size(ED_WeightMatrix))-ED_DegreeMatrix^(-1)*ED_WeightMatrix;


%% Compute Re-Scaled Euclidean Weight Matrix, ala Hoffmann

RescalingFactor=(ED_Degrees.*ED_Degrees').^(1-(q/2));
ED_WeightMatrix_Rescaled=ED_WeightMatrix./RescalingFactor;
ED_DegreeMatrix_Rescaled=diag(sum(ED_WeightMatrix_Rescaled,2));

L_ED_Rescaled=...
    ED_DegreeMatrix_Rescaled^((1-m)/(q-1))...
    *(ED_DegreeMatrix_Rescaled-ED_WeightMatrix_Rescaled)...
    *ED_DegreeMatrix_Rescaled^(-r/(q-1));

Time_Build_ED_Rescaled_L=toc;

%% Make sure the average degree weight is the same on the three graphs

%L_FD=trace(L_ED_Rescaled)/trace(L_FD)*L_FD;
%L_ED=trace(L_ED_Rescaled)/trace(L_ED)*L_ED;
%% Get spectra

tic;

[Vecs_FD, Vals_FD] = eigs(L_FD,NumEig,'smallestreal');
Vals_FD=diag(Vals_FD);

Time_FD_Eig=toc;

tic;

[Vecs_ED_Rescaled, Vals_ED_Rescaled] = eigs(L_ED_Rescaled,NumEig,'smallestreal');
Vals_ED_Rescaled=diag(Vals_ED_Rescaled);

Time_ED_Rescaled_Eig=toc;

[Vecs_ED, Vals_ED] = eigs(L_ED,NumEig,'smallestreal');
Vals_ED=diag(Vals_ED);

%% Compute Times

Time_FD{1}=Time_Euc_NN;
Time_FD{2}=Time_FD_NN;
Time_FD{3}=Time_Build_FD_L;
Time_FD{4}=Time_FD_Eig;

Time_ED_Rescaled{1}=Time_Euc_NN;
Time_ED_Rescaled{2}=Time_Build_ED_Rescaled_L;
Time_ED_Rescaled{3}=Time_ED_Rescaled_Eig;

%% Put Everything Into a Nice Struct

FD.Vals=Vals_FD;
FD.Vecs=Vecs_FD;
FD.Times=Time_FD;
FD.W=FD_WeightMatrix;
FD.D=FD_DegreeMatrix;
FD.L=L_FD;

ED_Rescaled.Vals=Vals_ED_Rescaled;
ED_Rescaled.Vecs=Vecs_ED_Rescaled;
ED_Rescaled.Times=Time_ED_Rescaled;
ED_Rescaled.W=ED_WeightMatrix_Rescaled;
ED_Rescaled.D=ED_DegreeMatrix_Rescaled;
ED_Rescaled.L=L_ED_Rescaled;

ED.Vals=Vals_ED;
ED.Vecs=Vecs_ED;
ED.W=ED_WeightMatrix;
ED.D=ED_DegreeMatrix;
ED.L=L_ED;


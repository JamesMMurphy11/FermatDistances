%Look at how different choices of p impact optimal graph cuts on a toy
%example that is highly elongated with a large density gap.

clear all;
close all;
profile off
profile on
SavePlots=0;

tic;

rng(42)

addpath(genpath('/Users/jmurph17/pwspd'));
setenv('PATH','/usr/bin:/usr/local/bin:/Library/TeX/texbin');

n=1000; 
Cratio=.05; 
beta=.1;
L=5;
X=[L,1].*(rand(n,2)-[.5,.5]);


BorderPoints=find(abs(X(:,2))<beta);
X(BorderPoints(rand(length(BorderPoints),1)<(1-Cratio)),:)=[];

GeoCut=ones(size(X,1),1);
GeoCut(X(:,1)>0)=2;
DensityCut=ones(size(X,1),1);
DensityCut(X(:,2)>0)=2;


%%

n=size(X,1);

figure;
scatter(X(:,1),X(:,2),[],GeoCut);
axis equal;
title('Geometry Cut','Interpreter','latex','FontSize',20);

figure;
scatter(X(:,1),X(:,2),[],DensityCut);
axis equal;
title('Density Cut','Interpreter','latex','FontSize',20);

%%

SigmaRange=[.01:.01:.1];

p_candidates=[1:.25:3]; % Weight for FD
NumEig=10;
NumClusters=2;

kNN=n;

[EucNN_Idx, EucNN_Distances] = knnsearch(X,X,'k',kNN);
Base=ones(n,kNN);

for j=1:n
    Base(j,:)=j*Base(j,:);
end

EucNN_Distances_T=EucNN_Distances';
EucNN_Idx_T=EucNN_Idx';
Base_T=Base';

A=sparse(Base_T(:),EucNN_Idx_T(:),EucNN_Distances_T(:),n,n);
ED = max(A, A');


for i=1:length(p_candidates)
    
    i

    p_FD=p_candidates(i);
    EucNN_Distances_p=EucNN_Distances.^p_FD;
    
    FD_NN_Idx = zeros(size(EucNN_Idx));
    FD_NN_Distances = zeros(size(EucNN_Distances));

    for k=1:n
        [temp_idx, temp_FD] = Dijkstra_with_early_stopping(EucNN_Idx, EucNN_Distances_p, k);
        FD_NN_Idx(k,:) = temp_idx';
        FD_NN_Distances(k,:) = temp_FD';
    end
    
    FD_NN_Distances_T=FD_NN_Distances';
    FD_NN_Idx_T=FD_NN_Idx';
    
    A_p=sparse(Base_T(:),FD_NN_Idx_T(:),FD_NN_Distances_T(:),n,n);
    FD{i}=max(A_p, A_p');

    for j=1:length(SigmaRange)
        
        sigma=SigmaRange(j);
                
        [L_FD_EigVecs{i,j},L_FD_EigVals{i,j},Labels_FD{i,j},K_FD(i,j),W_FD{i,j}, L_FD{i,j}, Time_FD(i,j)]=Embedding_PassW(FD{i},sigma,NumEig,NumClusters,'Gaussian');

        OA_Geo_FD(i,j) = GetAccuracies(Labels_FD{i,j},GeoCut,NumClusters);
        OA_Dens_FD(i,j) = GetAccuracies(Labels_FD{i,j},DensityCut,NumClusters);
       
        
        
    end
end
%%  Plot the two cuts

h=figure;
scatter(X(:,1),X(:,2),[],GeoCut);
axis equal;
title('Geometric Cut','Interpreter','latex','FontSize',20);

if SavePlots
    cd('/Users/jmurph17/JMLR2024_Code_Camera_Ready/JMLR2024_Code/DirichletEnergies/Images')
    saveas(h,'GeometryCut.pdf');
    system('pdfcrop --verbose GeometryCut.pdf');
    delete('GeometryCut.pdf');
end


h=figure;
scatter(X(:,1),X(:,2),[],DensityCut);
axis equal;
title('Density Cut','Interpreter','latex','FontSize',20);

if SavePlots
    cd('/Users/jmurph17/JMLR2024_Code_Camera_Ready/JMLR2024_Code/DirichletEnergies/Images')
    saveas(h,'DensityCut.pdf');
    system('pdfcrop --verbose DensityCut.pdf');
    delete('DensityCut.pdf');
end

%%  Plot SC accuracies w.r.t. different data partitions 

h=figure;
imagesc(OA_Geo_FD)
title('FD SC Accuracy w.r.t. Geometric Partition','Interpreter','latex','FontSize',20);
xlabel('$\epsilon$','Interpreter','latex','FontSize',20);
ylabel('$p$','Interpreter','latex','FontSize',20);
xticks(1:1:length(SigmaRange))
xticklabels(SigmaRange)
yticks(1:1:length(p_candidates));
yticks(1:1:length(p_candidates));
yticklabels(p_candidates);
colorbar

if SavePlots
    cd('/Users/jmurph17/JMLR2024_Code_Camera_Ready/JMLR2024_Code/DirichletEnergies/Images')
    saveas(h,'OA_Geo_FD.pdf');
    system('pdfcrop --verbose OA_Geo_FD.pdf');
    delete('OA_Geo_FD.pdf');
end

h=figure;
imagesc(OA_Dens_FD)
title('FD SC Accuracy w.r.t. Density Partition','Interpreter','latex','FontSize',20);
xlabel('$\epsilon$','Interpreter','latex','FontSize',20);
ylabel('$p$','Interpreter','latex','FontSize',20);
xticks(1:1:length(SigmaRange))
xticklabels(SigmaRange)
yticks(1:1:length(p_candidates));
yticklabels(p_candidates);
colorbar

if SavePlots
    cd('/Users/jmurph17/JMLR2024_Code_Camera_Ready/JMLR2024_Code/DirichletEnergies/Images')
    saveas(h,'OA_Dens_FD.pdf');
    system('pdfcrop --verbose OA_Dens_FD.pdf');
    delete('OA_Dens_FD.pdf');
end
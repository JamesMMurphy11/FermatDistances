% We conjecture the continuum operator that emerges from discrete Fermat
% distances matches on that emerges from a simpler construction
tic;

clearvars -except p thresh DataName
close all;

setenv('PATH','/usr/bin:/usr/local/bin:/Library/TeX/texbin');


profile off
profile on

SavePlots=1;

d=2;
NumTrials=10;
NumEig=100; 
N=5000;

%%  Make a folder for results

FolderName=['~/JMLR2024_Code/JMLR2024_Code/ComparisonExperiments/Results/',DataName,'/Example_p=',num2str(p),...
    '_thresh=',num2str(thresh),'_N=',num2str(N),'_T=',num2str(NumTrials)];

mkdir(FolderName)
cd(FolderName)
%%


for t=1:NumTrials

    if strcmp(DataName,'BallGap')
        X=SampleNonUniform_Dumbbell(N,thresh);
    elseif strcmp(DataName,'GaussianMixture')
        X=SampleGaussianMixture(N,thresh);
    end

    [FD{t}, ED_Rescaled{t}, ED{t}] = GetSpectra(X,p,NumEig);

    t

end

%%

for t=1:NumTrials
    Vals_FD(:,t)=FD{t}.Vals;
end

Vals_FD_Average=mean(Vals_FD,2);

for t=1:NumTrials
    Vals_ED_Rescaled(:,t)=ED_Rescaled{t}.Vals;
end

Vals_ED_Rescaled_Average=mean(Vals_ED_Rescaled,2);

for t=1:NumTrials
    Vals_ED(:,t)=ED{t}.Vals;
end

Vals_ED_Average=mean(Vals_ED,2);

%% Put the values on the same scale

Vals_FD_Average=sum(Vals_ED_Rescaled_Average(2:10))/sum(Vals_FD_Average(2:10)).*Vals_FD_Average;
Vals_ED_Average=sum(Vals_ED_Rescaled_Average(2:10))/sum(Vals_ED_Average(2:10)).*Vals_ED_Average;

%% Compute variances 

Vals_ED_Rescaled_Var=var(Vals_ED_Rescaled,0,2);
Vals_FD_Var=var(Vals_FD,0,2);
Vals_ED_Var=var(Vals_ED,0,2);



%% Plot eigenvalues

Width=2;

h=figure;
plot(Vals_FD_Average,'LineWidth',Width);
hold on;
plot(Vals_ED_Rescaled_Average,'LineWidth',Width);
hold on;
plot(Vals_ED_Average,'LineWidth',Width);
legend('Fermat Eigenvalues', 'Degree Renormalized Euclidean Eigenvalues','Euclidean Eigenvalues','Location','Northwest')
title(['Eigenvalues, $p$ = ',num2str(p), ', $\tau=$ ', num2str(thresh)],'Interpreter','latex','FontSize',18)

if SavePlots
    saveas(h,'Vals.pdf');
    system('pdfcrop --verbose Vals.pdf');
    delete('Vals.pdf');
end

% Plot eigenvalue standard deviations

h=figure; 
plot(sqrt(Vals_FD_Var(1:1:NumEig)),'LineWidth',Width);
hold on;
plot(sqrt(Vals_ED_Rescaled_Var(1:1:NumEig)),'LineWidth',Width);
hold on;
plot(sqrt(Vals_ED_Var(1:1:NumEig)),'LineWidth',Width);
axis tight
legend('Fermat Eigenvalues', 'Degree Renormalized Euclidean Eigenvalues','Euclidean Eigenvalues','Location','Southeast')
title(['Eigenvalues Standard Deviations, $p$ = ',num2str(p), ', $\tau=$ ', num2str(thresh)],'Interpreter','latex','FontSize',18)


if SavePlots
    saveas(h,'Vals_SD.pdf');
    system('pdfcrop --verbose Vals_SD.pdf');
    delete('Vals_SD.pdf');
end

% Plot eigenvalue absolute differences in log scale

h=figure;
plot(log10(abs(Vals_FD_Average(2:end)-Vals_ED_Rescaled_Average(2:end))),'LineWidth',Width);
title(['Log of Absolute Difference in Eigenvalues, $p$ = ',num2str(p), ', $\tau=$ ', num2str(thresh)],'Interpreter','latex','FontSize',18)

if SavePlots
    saveas(h,'Vals_LogDifferences.pdf');
    system('pdfcrop --verbose Vals_LogDifferences.pdf');
    delete('Vals_LogDifferences.pdf');
end

% Plot eigenvalue relative differences

h=figure;
plot((Vals_FD_Average(2:end)-Vals_ED_Rescaled_Average(2:end))./Vals_FD_Average(2:end),'LineWidth',Width);
title(['Relative Difference in Eigenvalues, $p$ = ',num2str(p), ', $\tau=$ ', num2str(thresh)],'Interpreter','latex','FontSize',18)

if SavePlots
    saveas(h,'Vals_RelDifferences.pdf');
    system('pdfcrop --verbose Vals_RelDifferences.pdf');
    delete('Vals_RelDifferences.pdf');
end

% Plot eigenvalue relative differences in log scale

h=figure;
plot(log10(abs((Vals_FD_Average(2:end)-Vals_ED_Rescaled_Average(2:end))./Vals_FD_Average(2:end))),'LineWidth',Width);
title(['Log of Relative Difference in Eigenvalues, $p$ = ',num2str(p), ', $\tau=$ ', num2str(thresh)],'Interpreter','latex','FontSize',18)


if SavePlots
    saveas(h,'Vals_LogRelDifferences.pdf');
    system('pdfcrop --verbose Vals_LogRelDifferences.pdf');
    delete('Vals_LogRelDifferences.pdf');
end


% Plot a representative data set

h=figure;
scatter(X(:,1),X(:,2))
title(['Data', ', $\tau=$ ', num2str(thresh), ', $n=$ ', num2str(N)],'Interpreter','latex','FontSize',18)


if SavePlots
    saveas(h,'Data.pdf');
    system('pdfcrop --verbose Data.pdf');
    delete('Data.pdf');
end

close all;

%%

save(['Results_p=',num2str(p),'_thresh=',num2str(thresh),'_N=',num2str(N),'_T=',num2str(NumTrials),'.mat'],'-v7.3')


%%

%{
figure;
plot(Vals_FD,'LineWidth',3);
hold on;
plot(Vals_ED,'LineWidth',3);
hold on;
plot(Vals_ED_Rescaled,'LineWidth',3);
legend('Fermat Eigenvalues', 'Euclidean Eigenvalues','Degree Renormalized Euclidean Eigenvalues')
%}


%legend('Fermat Eigenvalues', 'Degree Renormalized Euclidean Eigenvalues');

%{

%% Compute FD Laplacian with homogenous FD (i.e., adding 1/p power to exponent)

FD_Homo_WeightMatrix = sparse(Base_T(:),FD_IDX_T(:),(n)^(-1)*Sigma_FD_Homo^(-2).*exp(-(FD_KNN_Homo_T(:)./Sigma_FD_Homo).^2),n,n);
FD_Homo_WeightMatrix = max(FD_Homo_WeightMatrix, FD_Homo_WeightMatrix');
FD_Homo_WeightMatrix(logical(eye(size(FD_Homo_WeightMatrix))))=0;

FD_Homo_DegreeMatrix = diag(sum(FD_Homo_WeightMatrix));

L_FD_Homo = FD_Homo_DegreeMatrix - FD_Homo_WeightMatrix;

%figure;
%imagesc(FD_WeightMatrix)

[Vecs_FD_Homo, Vals_FD_Homo] = eigs(L_FD_Homo,NumEig,'smallestreal');

if Vecs_FD_Homo(1,2)<0
    Vecs_FD_Homo(:,2)=-Vecs_FD_Homo(:,2);
end

subplot(2,2,2);
scatter(X(:,1),X(:,2),[],Vecs_FD_Homo(:,2));
title('2nd Eigenfunction, FD Homo Laplacian');
colorbar


%% Plot difference in eigenvectors

EigVecDifferentLaplacians(1,:)=Vecs_ED(:,2);
EigVecDifferentLaplacians(2,:)=Vecs_ED_Rescaled(:,2);
EigVecDifferentLaplacians(3,:)=Vecs_FD(:,2);
EigVecDifferentLaplacians(4,:)=Vecs_FD_Homo(:,2);

DiffMatrix=squareform(pdist(EigVecDifferentLaplacians));

figure;
imagesc(DiffMatrix);
set(gca, 'XTick', [1:4], 'XTickLabel', {'Euc','Rescaled Euc','Fermat','Fermat Homo'}) % 10 ticks 
set(gca, 'YTick', [1:4], 'YTickLabel', {'Euc','Rescaled Euc','Fermat','Fermat Homo'}) % 10 ticks 
axis square;
colorbar

toc;

%%

if Vecs_ED(1,2)<0
    Vecs_ED(:,2)=-Vecs_ED(:,2);
end


subplot(2,2,4);
scatter(X(:,1),X(:,2),[],Vecs_ED(:,2));
title('2nd Eigenfunction, Euc. Laplacian');
colorbar

if Vecs_ED_Rescaled(1,2)<0
    Vecs_ED_Rescaled(:,2)=-Vecs_ED_Rescaled(:,2);
end

subplot(2,2,3);
scatter(X(:,1),X(:,2),[],Vecs_ED_Rescaled(:,2));
title('2nd Eigenfunction, Re-Scaled Euc. Laplacian');
colorbar

%%

%figure;
%imagesc(FD_WeightMatrix)


if Vecs_FD(1,2)<0
    Vecs_FD(:,2)=-Vecs_FD(:,2);
end
    
figure;
subplot(2,2,1);
scatter(X(:,1),X(:,2),[],Vecs_FD(:,2));
title('2nd Eigenfunction, FD Laplacian');
colorbar

[Vecs_ED, Vals_ED] = eigs(L_ED,NumEig,'smallestreal');

Sigma_FD_Homo=mean(FD_kNN_Homo(:));

% Plot eigenvalue differences

h=figure;
plot(abs(Vals_FD_Average-Vals_ED_Rescaled_Average),'LineWidth',Width);
title(['Absolute Difference in Eigenvalues, $p$ = ',num2str(p), ', $\tau=$ ', num2str(thresh)],'Interpreter','latex','FontSize',18)

if SavePlots
    saveas(h,'Vals_Differences.pdf');
    system('pdfcrop --verbose Vals_Differences.pdf');
    delete('Vals_Differences.pdf');
end

%}



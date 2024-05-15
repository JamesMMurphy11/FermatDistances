clear all;
close all;

p=2;  % density weight
s = 2; % normalization:

rgb =imread('Data/BlueSky/blue_sky_light_green_elongated_small.jpeg');
spatial_factor = 3;
patch_rad = 3;
sigma_KNN = 10; %nearest neighbor for sigma scaling

[rows, cols, ~] = size(rgb);
reduced_rgb = rgb(patch_rad+1:rows-patch_rad-1,patch_rad+1:cols-patch_rad-1,:);

 %% Make feature matrix

 [rrows, rcols, numberOfColorChannels] = size(reduced_rgb);
 X = zeros(rrows*rcols, numberOfColorChannels*(2*patch_rad+1)^2+2);
 for i=1:rrows
    for j=1:rcols
        rgb_row_idx = i+patch_rad;
        rgb_col_idx = j+patch_rad;
        sub_image = rgb( rgb_row_idx - patch_rad:rgb_row_idx + patch_rad, rgb_col_idx - patch_rad:rgb_col_idx + patch_rad,:);
        color_features = double(sub_image(:))/255;
        spatial_norm = sqrt(rrows^2+rcols^2);
        spatial_features = [rgb_row_idx rgb_col_idx]/spatial_norm;
        fc = length(color_features);
        X(rrows*(j-1)+i,:) =  [spatial_factor*(fc/(fc+2))*spatial_features (2/(fc+2))*color_features'];
     end
 end
 
 
 %% Make sparse kNN graph:
 
K=100;
n=size(X,1);

[IDX, D_KNN] = knnsearch(X,X,'k',K);

%% Compute FDs

tic
FD_IDX = zeros(size(IDX));
FD_KNN = zeros(size(D_KNN));
D_KNN_p = D_KNN.^p;
for i=1:n
    [temp_idx, temp_FD] = Dijkstra_with_early_stopping(IDX, D_KNN_p, i);
    temp_FD = temp_FD.^(1/p);
    FD_IDX(i,:) = temp_idx';
    FD_KNN(i,:) = temp_FD';
end
toc

%%
Base=ones(size(X,1),K);
for j=1:n
    Base(j,:)=j*Base(j,:);
end

D_KNN_T=D_KNN';
IDX_T=IDX';
FD_KNN_T=FD_KNN';
FD_IDX_T=FD_IDX';
Base_T=Base';

%% Compute Euclidean distance and Weights

sigma = mean(D_KNN(:,sigma_KNN));
Euc_dis = sparse(Base_T(:),IDX_T(:),D_KNN_T(:),n,n);
W = sparse(Base_T(:),IDX_T(:),exp(-D_KNN_T(:).^2/sigma^2),n,n);
Euc_dis = max(Euc_dis, Euc_dis');
W = max(W, W');

FD_sigma = mean(FD_KNN(:,sigma_KNN));

FD_dis = sparse(Base_T(:),FD_IDX_T(:),FD_KNN_T(:),n,n);
FD_W = sparse(Base_T(:),FD_IDX_T(:),exp(-FD_KNN_T(:).^2/FD_sigma^2),n,n);
FD_dis = max(FD_dis, FD_dis');
FD_W = max(FD_W, FD_W');

 %% Compute regular random walk Laplacian:
 
Deg = sparse(diag(sum(W,1)));
Deg_inv = sparse(diag(1./sum(W,1)));
tic
[V,D] = eigs(Deg-W,Deg,5,'smallestreal');
toc

 
 %% Compute FD Laplacian:
 
FD_Deg = sparse(diag(sum(FD_W,1)));
FD_W_weighted = sparse(size(FD_W));
if s==2
    FD_W_weighted = FD_W;
    FD_Deg_weighted = FD_Deg;
else
    alpha = 1-s/2;
    FD_Deg_alpha_inv = sparse(diag(1./(sum(FD_Deg,1).^alpha)));
    FD_W_weighted = FD_Deg_alpha_inv*FD_W*FD_Deg_alpha_inv;
    FD_Deg_weighted = sparse(diag(sum(FD_W_weighted,1)));
end
    
[FD_V,FD_D] = eigs(FD_Deg_weighted-FD_W_weighted,FD_Deg_weighted,5,'smallestreal');


%% Plot image

figure
imagesc(rgb);
axis equal
axis off

%% Plot just of regular v_2
figure
imagesc(reshape(V(:,2),size(reduced_rgb(:,:,1))))
axis equal
axis off


%% Plot of just Fermat v_2

figure
imagesc(reshape(FD_V(:,2),size(reduced_rgb(:,:,1))))
axis equal
axis off


 
 
clc; clear; close all
% Lecture 4 - PCA/MAF

addpath('..\Spectral Labs')
addpath('..\')

fid = fopen('avirisBand','r');
X_vec = fread(fid,'uint8');
fclose(fid);

nrows = 180;
ncols = 360;
nvars =  30;

% Binary file is loaded as one long vector, so we arrange it in 30 spectral
% images of 180 x 360 pixels. 
X = reshape(X_vec, ncols, nrows, nvars);
X = permute(X, [2 1 3]);

%% Plotting the 30 different bands
MakeFig(10, 10, 800, 400)
for i = 1:nvars
    subplot(5, 6, i)
    imagesc(X(:, :, i))
    colormap gray
    xlim([1, ncols])
    ylim([1, nrows])
    title(['Band nr. ',num2str(i),''])
end

%% Histograms
MakeFig(10, 10, 800, 400)
for i = 1:nvars
    subplot(5, 6, i)
    histogram(X(:, :, i), 50)
    title(['Hist for band nr. ',num2str(i),''])
end
return
%% PCA
% We reshape once again, so that we can do the statistical calculations
X_unmean = reshape(X, ncols*nrows, nvars);

% Column centering -> unmeaning
X_unmean = X_unmean - repmat(mean(X_unmean), ncols*nrows, 1);
% Variance covariance (dispersion) matrix
S = cov(X_unmean);
% Inverse of dispersion
iS = inv(S);

% Eigenvalues and eigenvectors
[V,Lambda] = eig(S);
lambda = diag(Lambda);

% Sorting eigen values in descendig order, and using index to sort eigen
% vectors likewise
[lambda, sort_idx] = sort(lambda,'descend');
V = V(:, sort_idx);

% Principal components
Q = X_unmean*V;
% Rearranging Principal Components
Q_im = reshape(Q, nrows, ncols, nvars);

MakeFig(10, 10, 800, 400)
for i = 1:nvars
    subplot(5, 6, i)
    imagesc(Q_im(:, :, i))
    colormap gray
    xlim([1, ncols])
    ylim([1, nrows])
    title(['PC nr. ',num2str(i),''])
end

% We check if our PCA is correct if cov(Q) = Lambda
isequal(round(diag(cov(Q)), 2), round(lambda, 2))

%% Checking variance proportions
% Proportion
prop = zeros(size(lambda));
for i = 1:length(prop)
    prop(i) = sum(lambda(1:i))/sum(lambda)*100;
end
format long g
disp([lambda,prop])
format short

%% Principal Components
% We see from above section that only the first to components are necessary
close all; clc

MakeFig(10, 10, 600, 200)
subplot(1, 2, 1)
imagesc(Q_im(:, :, 1))
colormap gray
axis equal
xlim([1, ncols])
ylim([1, nrows])
title(['PC nr. ',num2str(1),''])

subplot(1, 2, 2)
imagesc(Q_im(:, :, 2))
colormap gray
axis equal
xlim([1, ncols])
ylim([1, nrows])
title(['PC nr. ',num2str(2),''])
%% MAF
clc; close all
% We got Sigma already. Now we need compute Sigma_delta
% We use slicing for this
X_deltax = zeros(size(X));
X_deltay = zeros(size(X));

for i = 1:nvars
    % We compute the difference between a band, and the same band with all
    % the pixels shifted one position. For X_deltax we shift one to the
    % left, for the X_deltay we shift one up.
    
    % The band to subtract from
    band = X(:, :, i);
    
    % The deltax and deltay band 
    band_x = [band(:, 2:end), band(:, 1)];
    band_y = [band(2:end, :); band(1, :)];
    
    X_deltax(:, :, i) = band - band_x;
    X_deltay(:, :, i) = band - band_y;
end
% We reshape and column center the data
X_deltax = reshape(X_deltax, ncols*nrows, nvars);
X_deltay = reshape(X_deltay, ncols*nrows, nvars);
X_col = reshape(X, ncols*nrows, nvars);

% Column centering -> unmeaning
X_deltax = X_deltax - repmat(mean(X_deltax), ncols*nrows, 1);
X_deltay = X_deltay - repmat(mean(X_deltay), ncols*nrows, 1);

% % Normalizing data
% X_norm = normalize(X_col);
% X_deltax = normalize(X_deltax);
% X_deltay = normalize(X_deltay);


% Variance covariance (dispersion) matrix
S_norm = cov(X_unmean);
S_deltax = cov(X_deltax);
S_deltay = cov(X_deltay);

% Mean of the two dispersion matrices to get Disper
S_delta = 0.5*(S_deltax + S_deltay);

%% Checking our matrices
% Given that we have normalized the data, the dispersion matrix and
% correlation matrix should be the same
isequal(round(S_norm, 2), round(corr(X_col), 2))
% They are. Great


%% MAF - continued
% We will now find the Rayleigh Quotient
[U, Lambda] = eig(S_norm, S_delta);
lambda = diag(Lambda);
[lambda, sort_idx] = sort(lambda,'descend');
U = U(:, sort_idx);

% We get the projected eigenvalues as factors
F = X_unmean*U;
% Rearranging Factors
F_im = reshape(F, nrows, ncols, nvars);

%% Visualizing MAF
MakeFig(10, 10, 800, 400)
for i = 1:nvars
    subplot(5, 6, i)
    imagesc(F_im(:, :, i))
    colormap gray
    xlim([1, ncols])
    ylim([1, nrows])
    title(['Factor nr. ',num2str(i),''])
end

% Here we see how the first couple of factors tell us more about the data,
% whereas for the principal components, there some of the later components
% would still have less noise

% We check if our MAF is correcet by corr(F) = I
corr(F)
% Hooray, it is

%% One image
clc
BigIm = zeros(6*nrows, 5*ncols);
idx_count = 0;
for i = 1:6
    for j = 1:5
        idx_count = idx_count + 1;
        BigIm((1 + (i-1)*nrows):i*nrows, (1 + (j-1)*ncols):j*ncols) = F_im(:, :, idx_count);
    end
end

MakeFig(10, 10, 800, 400)
imagesc(BigIm)
colormap gray
axis equal
xlim([1, 5*ncols])
ylim([1, 6*nrows])
%% 
MakeFig(10, 10, 800, 400)
% Standard deviations in which to show the respective image

imshowrgb(X, [2, 4, 6], 2)
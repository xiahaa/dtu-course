clc; clear; close all
% Lecture 5 - Supervised Classification SVM

% Adds a path into our data folder
addpath('..\Spectral Labs')
addpath('..\') % For the MakeFig.m function

fid = fopen('igalmss','r');
X_vec = uint8(fread(fid,'uint8').*2'^2);
fclose(fid);

% From igalmss.hdr we get dimensions
nrows = 512;
ncols = 512;
nvars = 4;
nobs = ncols*nrows;

% We arrange the matrix so we have each pixels values in all spectra
% side-by-side
% X as a matrix
X_mat = double(reshape(X_vec, nobs, nvars));
% For displaying image
X = reshape(X_vec, ncols, nrows, nvars);

imshowrgb(X,[4, 2, 1])

%% PCA
% Since the SVM is very computationally heavy, we will perform PCA in order
% to reduce the amounts of data input

% Column centering -> unmeaning
X_unmean = X_mat - repmat(mean(X_mat), ncols*nrows, 1);
% Variance covariance (dispersion) matrix
S = cov(X_unmean);

% Eigenvalues and eigenvectors
[V, Lambda] = eig(S);
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
    subplot(2, 2, i)
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

% We see that including just the first 2 principal components will give us
% all more than 98% of the variance in the data
%% Training images
% close all; clc
% % Reading of interest polygon
% t1 = roipoly(X(:,:,1:3)); % water
% t2 = 2*roipoly(X(:,:,1:3)); % everything else
% 
% X_train = t1+t2;
% clear t1 t2
% close all

%% Random sample
% Load already determined training areas
load training_SVM.mat
close all; clc
% We now got our traning image from which we can extract to data sets of
% pixels corresponding to the defined classes
data_1 = Q(X_train==1,:);
data_2 = Q(X_train==2,:);

% Number of pixels in each training image
N_1 = size(data_1, 1);
N_2 = size(data_2, 1);

rng(1) % For reproducibility

% We collect random samples
N_samp = 1000; % Number of samples
% First we get indices
idx_1 = randsample(1:N_1, N_samp);
idx_2 = randsample(1:N_2, N_samp);
% Then samples
samp_1 = double(data_1(idx_1, :));
samp_2 = double(data_2(idx_2, :));

samples = [samp_1; samp_2];
samples = samples(:, 1:2);

% Plotting the samples
MakeFig(100, 50, 800, 400)
plot(samp_1(:, 1), samp_1(:, 2), 'r.','MarkerSize',15)
hold on
plot(samp_2(:, 1), samp_2(:, 2), 'b.','MarkerSize',15)
axis equal
xlabel('PC 1')
ylabel('PC 2')
legend('Sample 1', 'Sample 2')
% We see that most of the points are clearly distinguished. The water is
% seemingly clearly defined, except for a few outliers among the 1000
% samples. Everything else is of course a bit more vague.

% We apply class labels
classes = ones(N_samp*2, 1);
classes(1:N_samp) = -1;
return

%% SVM using Sigmoid Kernel
close all; clc
% Try adjusting gamma parameter in "mysigmoid.m". We find that gamma = 0.5
% is a good value
Mdl_1 = fitcsvm(samples, classes,'KernelFunction','mysigmoid','Standardize',true);

% Predict scores over the grid
d = 0.5; % Smaller grid will improve missclassification rate, but increase computation
[x1_Grid, x2_Grid] = meshgrid(min(samples(:,1)):d:max(samples(:,1)),min(samples(:,2)):d:max(samples(:,2)));
x_Grid = [x1_Grid(:),x2_Grid(:)];
[~,scores_1] = predict(Mdl_1,x_Grid); % The scores

MakeFig(100, 50, 800, 400)
h(1:2) = gscatter(samples(:,1), samples(:,2), classes);
hold on
h(3) = plot(samples(Mdl_1.IsSupportVector,1),...
    samples(Mdl_1.IsSupportVector,2),'ko','MarkerSize',10);
    % Support vectors
contour(x1_Grid,x2_Grid,reshape(scores_1(:,2),size(x1_Grid)),[0 0],'k');
    % Decision boundary
title('Scatter Diagram with the Decision Boundary')
legend({'-1','1','Support Vectors'},'Location','Best');
hold off

% Missclassification rate
miss_class = kfoldLoss(crossval(Mdl_1))
disp(['Missclassification rate is ',num2str(100*miss_class),' %'])

%% Predicting classes of remaining pixels
clc
[class_1, score_1] = predict(Mdl_1, Q(:, 1:2));

%% Visualizing
MakeFig(100, 50, 800, 400)
subplot(1, 2, 1)
imshowrgb(X,[1, 2, 3])
title('Original image displayed using band 1, 2 and 3')

subplot(1, 2, 2)
imshow(reshape(class_1, ncols, nrows))
colormap(jet)
title(' Classification using kernel')

%% Linear SVM
% Try adjusting gamma parameter in "mysigmoid.m". We find that gamma = 0.5
% is a good value
Mdl_2 = fitcsvm(samples, classes,'KernelFunction','rbf');

% Predict scores over the grid
d = 0.5; % Smaller grid will improve missclassification rate, but increase computation
[x1_Grid, x2_Grid] = meshgrid(min(samples(:,1)):d:max(samples(:,1)),min(samples(:,2)):d:max(samples(:,2)));
x_Grid = [x1_Grid(:),x2_Grid(:)];
[~,scores_1] = predict(Mdl_2,x_Grid); % The scores

MakeFig(100, 50, 800, 400)
h(1:2) = gscatter(samples(:,1), samples(:,2), classes);
hold on
h(3) = plot(samples(Mdl_2.IsSupportVector,1),...
    samples(Mdl_2.IsSupportVector,2),'ko','MarkerSize',10);
    % Support vectors
contour(x1_Grid,x2_Grid,reshape(scores_1(:,2),size(x1_Grid)),[0 0],'k');
    % Decision boundary
title('Scatter Diagram with the Decision Boundary')
legend({'-1','1','Support Vectors'},'Location','Best');
hold off

% Missclassification rate
miss_class = kfoldLoss(crossval(Mdl_2))
disp(['Missclassification rate is ',num2str(100*miss_class),' %'])

%% Predicting and Visualizing
clc
[class_2, score_2] = predict(Mdl_2, Q(:, 1:2));

MakeFig(100, 50, 800, 400)
subplot(1, 2, 1)
imshowrgb(X,[1, 2, 3])
title('Original image displayed using band 1, 2 and 3')

subplot(1, 2, 2)
imshow(reshape(class_2, ncols, nrows))
colormap(jet)
title('Classification with linear SVM')
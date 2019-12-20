clc;close all;clear all;

addpath ../../utils;

baseDir = '../../data/optical_flow_data/Basketball';

% load database
buildingScene = imageDatastore(baseDir);
numImages = numel(buildingScene.Files);

% load a certain image
im1 = imread(buildingScene.Files{1});
im1 = imPreprocessing(im1);

im2 = imread(buildingScene.Files{2});
im2 = imPreprocessing(im2);

% image size
% [flow_u, flow_v] = HS(im1, im2, 0.01, 1000);
[flow_u, flow_v] = HornSchunck(im1, im2, 1, 100);


% display dense flow as an image
img = computeColor(flow_u,flow_v);
figure
imshow(img, 'InitialMagnification',1000);hold on;

function showFlowQuiver(im, flow_u, flow_v)
    figure;imshow(im, 'InitialMagnification',1000);hold on;
    height = size(im,1);
    width = size(im,2);
    % opflow = opticalFlow(flow_u,flow_v);
    % plot(opflow,'DecimationFactor',[1 1],'ScaleFactor',1);
    [yy,xx] = meshgrid(1:height,1:width);
    xx = vec(xx');
    yy = vec(yy');
    quiver(xx,yy,flow_u(:),flow_v(:),'LineWidth',1.5, 'Color','r','MaxHeadSize',1);axis image
end

function [Ex, Ey, Et] = derivative(Im1, Im2)
    % function DERIVATIVE computes partial derivatives Ex, Ey, Et
    % of a sequence of 2 images Im1, Im2 of double class.

%     kernels for convolution
    Kx = 0.25 * [-1 1; -1 1];
    Ky = 0.25 * [-1 -1; 1 1];
    Kt = 0.25 * [-1 -1; -1 -1]; % kt1 = Kt, kt2 = -Kt

    % compute derivatives
    Ex = conv2(Im1, Kx, 'same') + conv2(Im2, Kx, 'same');
    Ey = conv2(Im1, Ky, 'same') + conv2(Im2, Ky, 'same');
    Et = conv2(Im1, Kt, 'same') + conv2(Im2, -Kt, 'same');

%     x = -1:1:1;
%     cons1 = 1*1;
%     hg = 1/sqrt(2*pi*cons1).*exp(-x.^2./(2*cons1)); hg = hg ./ sum(abs(hg(:)));
%     hgx = 1/sqrt(2*pi*cons1).*exp(-x.^2./(2*cons1)).*(-x./cons1); hgx = hgx ./ sum(abs(hgx(:)));
%     hgx = fliplr(hgx);% convolution kernel is the flip version.
%     Ex = imfilter(Im1,hgx,'replicate','same');
%     Ex = imfilter(Ex,hg','replicate','same');
%     Ey = imfilter(Im1,hgx','replicate','same');
%     Ey = imfilter(Ey,hg,'replicate','same');
%     Et = Im2 - Im1;
end

function smoothedImg=smoothImg(img,segma)
% Convolving an image with a Gaussian kernel.

% The degree of smoothing is determined by the Gaussian's standard
% deviation 'segma'. The Gaussian outputs a 'weighted average' of each
% pixel's neighborhood, with the average weighted more towards the value of
% the central pixels. The larger the standard deviation, the less weight
% towards the central pixels and more weight towards the pixels away, hence
% heavier blurring effect and less details in the output image.
%
% Author: Mohd Kharbat at Cranfield Defence and Security
% mkharbat(at)ieee(dot)org , http://mohd.kharbat.com
% Published under a Creative Commons Attribution-Non-Commercial-Share Alike
% 3.0 Unported Licence http://creativecommons.org/licenses/by-nc-sa/3.0/
%
% October 2008

if nargin<2
    segma=1;
end

G=gaussFilter(segma);
smoothedImg=conv2(img,G,'same');
smoothedImg=conv2(smoothedImg,G','same');
end


function G=gaussFilter(segma,kSize)
% Creates a 1-D Gaussian kernel of a standard deviation 'segma' and a size
% of 'kSize'. 
%
% In theory, the Gaussian distribution is non-zero everywhere. In practice,
% it's effectively zero at places further away from about three standard
% deviations. Hence the reason why the kernel is suggested to be truncated
% at that point.
%
% The 2D Gaussian filter is a complete circular symmetric operator. It can be
% seperated into x and y components. The 2D convolution can be performed by
% first convolving with 1D Gaussian in the x direction and the same in the
% y direction.
%
% Author: Mohd Kharbat at Cranfield Defence and Security
% mkharbat(at)ieee(dot)org , http://mohd.kharbat.com
% Published under a Creative Commons Attribution-Non-Commercial-Share Alike
% 3.0 Unported Licence http://creativecommons.org/licenses/by-nc-sa/3.0/
%
% October 2008

if nargin<1
    segma=1;
end
if nargin<2
    kSize=2*(segma*3);
end

x=-(kSize/2):(1+1/kSize):(kSize/2);
G=(1/(sqrt(2*pi)*segma)) * exp (-(x.^2)/(2*segma^2));
end

function [fx, fy, ft] = computeDerivatives(im1, im2)

% derivatives as in Barron
% fx= conv2(im1,(1/12)*[-1 8 0 -8 1],'same');
% fy= conv2(im1,(1/12)*[-1 8 0 -8 1]','same');
% ft = conv2(im1, 0.25*ones(2),'same') + conv2(im2, -0.25*ones(2),'same');
% fx=-fx;fy=-fy;

% An alternative way to compute the spatiotemporal derivatives is to use simple finite difference masks.
% fx = conv2(im1,[1 -1]);
% fy = conv2(im1,[1; -1]);
% ft= im2-im1;

if size(im2,1)==0
    im2=zeros(size(im1));
end

Xkern = 0.25* [-1 1; 
               -1 1];
               
Ykern = 0.25*[-1 -1; 
               1 1];
               
Tkern = 0.25*ones(2);

% Horn-Schunck original method
fx = conv2(im1,Xkern,'same') + conv2(im2, Xkern,'same');
fy = conv2(im1, Ykern, 'same') + conv2(im2, Ykern, 'same');
ft = conv2(im1, Tkern,'same') + conv2(im2, -Tkern, 'same');

end

function plotFlow(u, v, imgOriginal, rSize, scale)
% Creates a quiver plot that displays the optical flow vectors on the
% original first frame (if provided). See the MATLAB Function Reference for
% "quiver" for more info.
%
% Usage:
% plotFlow(u, v, imgOriginal, rSize, scale)
%
% u and v are the horizontal and vertical optical flow vectors,
% respectively. imgOriginal, if supplied, is the first frame on which the
% flow vectors would be plotted. use an empty matrix '[]' for no image.
% rSize is the size of the region in which one vector is visible. scale
% over-rules the auto scaling.
%
% Author: Mohd Kharbat at Cranfield Defence and Security
% mkharbat(at)ieee(dot)org , http://mohd.kharbat.com
% Published under a Creative Commons Attribution-Non-Commercial-Share Alike
% 3.0 Unported Licence http://creativecommons.org/licenses/by-nc-sa/3.0/
%
% October 2008
% Rev: Jan 2009

figure();

if nargin>2
    imshow(imgOriginal);
    hold on;
end
if nargin<4
    rSize=5;
end
if nargin<5
    scale=3;
end

% Enhance the quiver plot visually by showing one vector per region
for i=1:size(u,1)
    for j=1:size(u,2)
        if floor(i/rSize)~=i/rSize || floor(j/rSize)~=j/rSize
            u(i,j)=0;
            v(i,j)=0;
        end
    end
end
quiver(u, v, scale)%, 'color', 'b', 'linewidth', 2);
%set(gca,'YDir','reverse');
end
    
function [u, v] = HornSchunck(im1, im2, alpha, ite, uInitial, vInitial, displayFlow, displayImg)
% Horn-Schunck optical flow method 
% Horn, B.K.P., and Schunck, B.G., Determining Optical Flow, AI(17), No.
% 1-3, August 1981, pp. 185-203 http://dspace.mit.edu/handle/1721.1/6337
%
% Usage:
% [u, v] = HS(im1, im2, alpha, ite, uInitial, vInitial, displayFlow)
% For an example, run this file from the menu Debug->Run or press (F5)
%
% -im1,im2 : two subsequent frames or images.
% -alpha : a parameter that reflects the influence of the smoothness term.
% -ite : number of iterations.
% -uInitial, vInitial : initial values for the flow. If available, the
% flow would converge faster and hence would need less iterations ; default is zero. 
% -displayFlow : 1 for display, 0 for no display ; default is 1.
% -displayImg : specify the image on which the flow would appear ( use an
% empty matrix "[]" for no image. )
%
% Author: Mohd Kharbat at Cranfield Defence and Security
% mkharbat(at)ieee(dot)org , http://mohd.kharbat.com
% Published under a Creative Commons Attribution-Non-Commercial-Share Alike
% 3.0 Unported Licence http://creativecommons.org/licenses/by-nc-sa/3.0/
%
% October 2008
% Rev: Jan 2009

%% Default parameters
if nargin<1 || nargin<2
    im1=imread('yos9.tif');
    im2=imread('yos10.tif');
end
if nargin<3
    alpha=1;
end
if nargin<4
    ite=100;
end
if nargin<5 || nargin<6
    uInitial = zeros(size(im1(:,:,1)));
    vInitial = zeros(size(im2(:,:,1)));
elseif size(uInitial,1) ==0 || size(vInitial,1)==0
    uInitial = zeros(size(im1(:,:,1)));
    vInitial = zeros(size(im2(:,:,1)));
end
if nargin<7
    displayFlow=1;
end
if nargin<8
    displayImg=im1;
end
 
%% Convert images to grayscale
if size(size(im1),2)==3
    im1=rgb2gray(im1);
end
if size(size(im2),2)==3
    im2=rgb2gray(im2);
end
im1=double(im1);
im2=double(im2);

im1=smoothImg(im1,1);
im2=smoothImg(im2,1);

tic;

%%
% Set initial value for the flow vectors
u = uInitial;
v = vInitial;

% Estimate spatiotemporal derivatives
[fx, fy, ft] = computeDerivatives(im1, im2);

% Averaging kernel
kernel_1=[1/12 1/6 1/12;1/6 0 1/6;1/12 1/6 1/12];

% Iterations
for i=1:ite
    % Compute local averages of the flow vectors
    uAvg=conv2(u,kernel_1,'same');
    vAvg=conv2(v,kernel_1,'same');
    % Compute flow vectors constrained by its local average and the optical flow constraints
    der = ( fx.*uAvg + fy.*vAvg + ft ) ./ ( alpha^2 + fx.^2 + fy.^2);
    u = uAvg - fx .* der; 
    v = vAvg - fy .* der;
end

u(isnan(u))=0;
v(isnan(v))=0;

%% Plotting
if displayFlow==1
    plotFlow(u, v, displayImg, 5, 5); 
end

showFlowQuiver(im1, u, v);
end

function [U, V] = HS(Im1, Im2, alpha, N)
    % function HS computes flow velocities U, V of a sequence of 2 images
    % Im1, Im2 of double class based on Horn-Schunck algorithm. alpha is the
    % weighting factor and N is the number of iteration.

    % compute partial derivatives
    [Ex, Ey, Et] = derivative(Im1, Im2);

    % intial U, V
    [l, c] = size(Im1);
    U = zeros(l, c);
    V = zeros(l, c);

    K = [1/12 1/6 1/12; 1/6 0 1/6; 1/12 1/6 1/12]; % Laplacian kernel
    A = alpha^2 + Ex.^2 + Ey.^2;

    for i = 1:N
        % compute U,V averages
        U_avg = conv2(U, K, 'same');
        V_avg = conv2(V, K, 'same');
        B = (Ex.*U_avg + Ey.*V_avg + Et);

        % compute U, V at current iteration
        U = U_avg - Ex.*B./A;
        V = V_avg - Ey.*B./A;
    end
    showFlowQuiver(Im1, U, V);
end

function [flow_u, flow_v] = LK(Im1, Im2, ws)
    % function LK computes flow velocities U, V of a sequence of 2 images
    % Im1, Im2 of double class based on Lucas-Kanade method. ws is the
    % window size.
    % return (U, V) flow velocities, C coordinates of corners in Im1

    % dense flow
    width = size(Im1,2);
    height = size(Im1,1);
    [yy,xx] = meshgrid(1:height,1:width);
    xx = vec(xx');
    yy = vec(yy');
    C = [xx yy];

    % mark corners in the margin as (-1, -1)
    [m, n] = size(Im1);
    w = round(ws/2); % half of windows size
    for i = 1:size(C, 1)
        x = C(i, 1);
        y = C(i, 2);
        if (x <= w) || (x >= n-w) % corners in left, right margin
            C(i,:) = -1;
        end
        if (y <= w) || (y >= m-w) % corner in top, bottom margin
            C(i,:) = -1;
        end
    end

    % remove corners in margins
    id = C(:,1) == -1 | C(:,2) == -1;
    C(id,:) = [];

    % Lucas-Kanade constraint
    nc = size(C, 1); % number of suitable corners
    U = zeros(nc, 1);
    V = zeros(nc, 1);

    % compute partial derivatives
    [Ex, Ey, Et] = derivative(Im1, Im2);

    for i = 1:nc
        % get derivatives of neighborhood points
        x = C(i, 1); y = C(i, 2);
        Ix = Ex(y-w:y+w, x-w:x+w)';
        Iy = Ey(y-w:y+w, x-w:x+w)';
        It = Et(y-w:y+w, x-w:x+w)';

        % get A, b
        Ix = Ix(:); Iy = Iy(:); It = It(:);
        A = [Ix Iy];
        b = -It;

        % compute U, V
        X = pinv(A) * b;
        U(i) = X(1);
        V(i) = X(2);
    end
    flow_u = zeros(m,n);
    flow_v = zeros(m,n);
    flow_u(~id) = U;
    flow_v(~id) = V;
    showFlowQuiver(Im1, flow_u, flow_v);
end


function [flow_u, flow_v] = denseflowLK(im1, im2, iu, iv, hsize)
% compute optical flow using Lucas-Kanade
    % grad    
    [Ix, Iy] = grad2(im1,1,1);
    
    % dense flow
    width = size(im1,2);
    height = size(im1,1);
    [yy,xx] = meshgrid(1:height,1:width);
    xx = vec(xx');
    yy = vec(yy');
    
    if isempty(iu) 
        flow_u = zeros(height, width);
    else
        flow_u = iu;
    end
    
    if isempty(iv) 
    	flow_v = zeros(height, width);
    else
        flow_v = iv;
    end
    
    % precomputing
    Ix2 = Ix.*Ix;
    Iy2 = Iy.*Iy;
    Ixy = Ix.*Iy;
    % filtering
    ker = gauker(2);
    
    sIx2 = imfilter(Ix2,ker,'replicate','same');
    sIy2 = imfilter(Iy2,ker,'replicate','same');
    sIxy = imfilter(Ixy,ker,'replicate','same');
    sIxy2 = sIxy.*sIxy;

    % interpolation
%     It = gradIntensity(xx,yy,vec(flow_u),vec(flow_v),im1,im2);
%     It = reshape(It, height, width);
    It = im2 - im1;
    
    Ixt = Ix.*It;
    Iyt = Iy.*It;
        
    sIxt = imfilter(Ixt,ker,'replicate','same');
    sIyt = imfilter(Iyt,ker,'replicate','same');
    
    % LK-flow    
    % opt1: vectorization, faster
    % conditioning: check if the smallest eigenvalue is very close to zero.
    % since it is a 2x2 positive semidefinite matrix, its eigen value has 
    % a analytical solution and must be a real value greater or equal to 0.
%     eig_smallest = 0.5.*(sIx2 + sIy2 - sqrt((sIx2-sIy2).^2+4.*sIxy2));
%     eig_largest = 0.5.*(sIx2 + sIy2 + sqrt((sIx2-sIy2).^2+4.*sIxy2));
%     ratio = eig_smallest ./ eig_largest;
%     invalid = ratio < 0.01 | (eig_smallest < 1e-6);
    
    % second option is to use Harris 
%     score = sIx2.*sIy2 - sIxy2 - 0.01.*(sIx2+sIy2).^2;
    
    % third option is to use the hormonic mean
    score = (sIx2.*sIy2 - sIxy2)./(sIx2+sIy2);
    invalid = score < 0.1*max(score(:));
    
    % add a smaller value to the diagonal elements of those invalid pixels.
    %     sIx2(invalid) = sIx2(invalid) + 0.1;
    %     sIy2(invalid) = sIy2(invalid) + 0.1;

    % analytical inversion
    s = 1./(sIx2.*sIy2 - sIxy2);
    dflow_u = -( sIy2.*sIxt - sIxy.*sIyt).*s;
    dflow_v = -(-sIxy.*sIxt + sIx2.*sIyt).*s;
        
    dflow_u(invalid) = 0;
    dflow_v(invalid) = 0;
    
    flow_u = flow_u + dflow_u;
    flow_v = flow_v + dflow_v;
    
    % debug only
    showFlowQuiver(im1, flow_u, flow_v);
end

function [flow_u, flow_v] = denseflowHS(im1, im2, hsize)
% compute optical flow using Horn-Shunck method
   [flow_u, flow_v] = denseflowPyrLK(im1, im2);

    hsize = 3;
    sigma = 1;
    %-------- this ends the initialization, then start HS iteration --------%
    
    % grad    
    [Ix, Iy] = grad2(im1,3,1);
    
    % dense flow
    width = size(im1,2);
    height = size(im1,1);
    [yy,xx] = meshgrid(1:height,1:width);
    xx = vec(xx');
    yy = vec(yy');
    
    if isempty(iu) 
        flow_u = zeros(height, width);
    else
        flow_u = iu;
    end
    
    if isempty(iv) 
    	flow_v = zeros(height, width);
    else
        flow_v = iv;
    end
    
    % precomputing
    Ix2 = Ix.*Ix;
    Iy2 = Iy.*Iy;
    Ixy = Ix.*Iy;    
    
    % kernel 1, hard code 3x3 averaging
%     ker_avg = [1/12 1/6 1/12;1/6 0 1/6;1/12 1/6 1/12];
    ker_avg = gaussian_kernel_calculator(2, hsize, sigma);
    
    max_iter = 50;
    alpha = 0.5;
    for iter = 1:max_iter
        % arveraging
        ubar = imfilter(flow_u,ker_avg,'replicate','same');
        vbar = imfilter(flow_v,ker_avg,'replicate','same');
        % update
        den = alpha*alpha + Ix2 + Iy2;
        flow_u = ubar - (Ix2.*ubar + Ixy.*vbar + Ixt)./den;
        flow_v = vbar - (Ixy.*ubar + Iy2.*vbar + Iyt)./den;
    end
    
end


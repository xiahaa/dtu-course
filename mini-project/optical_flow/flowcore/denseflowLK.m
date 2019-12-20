function [flow_u, flow_v] = denseflowLK(im1, im2, iu, iv, verbose)
% compute optical flow using Lucas-Kanade
    % dense flow
    width = size(im2,2);
    height = size(im2,1);
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
    
    im2_warp = warpImage(im1, flow_u, flow_v, im2);
    % grad 
    [Ix, Iy] = grad2(im2_warp);
    It = im2_warp - im1;
    
    % precomputing
    Ix2 = Ix.*Ix;
    Iy2 = Iy.*Iy;
    Ixy = Ix.*Iy;
    % filtering
    ker = gauker(2, 1);% hsize + sigma
    
    sIx2 = imfilter(Ix2,ker,'symmetric','same');
    sIy2 = imfilter(Iy2,ker,'symmetric','same');
    sIxy = imfilter(Ixy,ker,'symmetric','same');
    sIxy2 = sIxy.*sIxy;
    
    Ixt = Ix.*It;
    Iyt = Iy.*It;

    sIxt = imfilter(Ixt,ker,'symmetric','same');
    sIyt = imfilter(Iyt,ker,'symmetric','same');
    
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
%     score = sIx2.*sIy2 - sIxy2 - 0.05.*(sIx2+sIy2).^2;
    
    % third option is to use the hormonic mean
    score = (sIx2.*sIy2 - sIxy2)./(sIx2+sIy2);
    det1 = sIx2.*sIy2 - sIxy2;
    invalid =  score < 0.05*max(score(:)) | abs(det1) < 1e-6;%
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
    
    if verbose == true
        % debug only
        showFlowQuiver(im1, flow_u, flow_v);
    end
end
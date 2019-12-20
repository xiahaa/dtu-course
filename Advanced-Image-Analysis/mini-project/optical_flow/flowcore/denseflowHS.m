
function [flow_u, flow_v] = denseflowHS(im1, im2, iu, iv, verbose)
% compute optical flow using Horn-Shunck method  
    height = size(im2,1);
    width = size(im2,2);
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
	
    % kernel 1, hard code 3x3 averaging
    %     ker_avg = [1/12 1/6 1/12;1/6 0 1/6;1/12 1/6 1/12];
	sigma = 1;
    hsize = 2;
    ker_avg = gauker(hsize, sigma);% hsize + sigma

    maxiter = 3;
    for iter = 1:maxiter
        im2_warp = warpImage(im1, flow_u, flow_v, im2);
        % grad    
        [Ix, Iy] = grad2(im2_warp);
        It = im2_warp - im1;
    
        % precomputing
        Ix2 = Ix.*Ix;
        Iy2 = Iy.*Iy;
        Ixy = Ix.*Iy;    

   	    sIx2 = imfilter(Ix2,ker_avg,'symmetric','same');
        sIy2 = imfilter(Iy2,ker_avg,'symmetric','same');
        sIxy = imfilter(Ixy,ker_avg,'symmetric','same');
        Ixt = Ix.*It;
        Iyt = Iy.*It;
        sIxt = imfilter(Ixt,ker_avg,'symmetric','same');
        sIyt = imfilter(Iyt,ker_avg,'symmetric','same');
        
        max_iter = 10;
        alpha = 1;
        
        uu = zeros(size(im2));
        vv = zeros(size(im2));
        
        for i = 1:max_iter
            % arveraging
            ubar = imfilter(uu,ker_avg,'replicate','same');
            vbar = imfilter(vv,ker_avg,'replicate','same');
            % update
            den = alpha*alpha + sIx2 + sIy2;
            
            du = (sIx2.*ubar + sIxy.*vbar + sIxt)./den;
            dv = (sIxy.*ubar + sIy2.*vbar + sIyt)./den;
            
            uu = ubar - du;
            vv = vbar - dv;
            
            if max(abs(du(:))) < 1e-3 && max(abs(dv(:))) < 1e-3
                disp('Exit 1~~~~~~');
                break;
            end
        end
        flow_u = flow_u + uu;
        flow_v = flow_v + vv;
        if max(abs(uu(:))) < 1e-3 && max(abs(vv(:))) < 1e-3
            disp('Exit 2~~~~~~');
            break;
        end
    end
    if verbose == true
        % debug only
        showFlowQuiver(im1, flow_u, flow_v);
    end
end
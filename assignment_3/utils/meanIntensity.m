function [cin, cout] = meanIntensity(im, curve, boundary)
%Compute mean intensity for deformable models.
    polygon = curve;
    polygon(:,end+1) = curve(:,1);
    % region inside curve
    mask = poly2mask(polygon(2,:),polygon(1,:),size(im,1),size(im,2));
    % if boundary is empty consider all pixels inside region as in, otherwise as out
    if isempty(boundary)
        maskin = mask;
        maskout = ~mask;
    else
        % otherwise, find relevent in and out
        maskin = mask & boundary;
        maskout = xor(maskin, boundary);
    end
    
    % debug, visualization
%     im(mask) = 1;
%     imshow(im);
    cin = mean(vec(im(maskin)));
    cout = mean(vec(im(maskout)));
end
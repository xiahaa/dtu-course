function [cin, cout] = meanIntensity(im, curve, boundary)
    polygon = curve;
    polygon(:,end+1) = curve(:,1);
    mask = poly2mask(polygon(2,:),polygon(1,:),size(im,1),size(im,2));
    
    if isempty(boundary)
        maskin = mask;
        maskout = ~mask;
    else
        maskin = mask & boundary;
        maskout = ~mask & boundary;
    end
    
    % debug
%     im(mask) = 1;
%     imshow(im);
    cin = mean(vec(im(maskin)));
    cout = mean(vec(im(maskout)));
end
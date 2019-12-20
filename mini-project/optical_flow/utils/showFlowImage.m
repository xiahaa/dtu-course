function showFlowImage(flow_u, flow_v)
    % display dense flow as an image
    img = computeColor(flow_u,flow_v);
    figure
    if size(img,1)*size(img,2) < 10000
        imshow(img, 'InitialMagnification',10000);hold on;
    else
        imshow(img);
    end
end

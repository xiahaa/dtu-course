function showFlowQuiver(im, flow_u, flow_v)
    figure;
    if size(im,1)*size(im,2) < 10000
        imshow(im, 'InitialMagnification',1000);hold on;
    else
        imshow(im);hold on;
    end
    height = size(im,1);
    width = size(im,2);
    
    [xx1,yy1] = meshgrid(1:width,1:height);
    
    xxs = linspace(1,width, 50);
    yys = linspace(1,height, 50);
    
    [xx2,yy2] = meshgrid(xxs, yys);
    uu = interp2(xx1,yy1,flow_u, xx2, yy2);
    vv = interp2(xx1,yy1,flow_v, xx2, yy2);
    
%     opflow = opticalFlow(uu,vv);
%     plot(opflow,'DecimationFactor',[1 1],'ScaleFactor',1);
%     [yy,xx] = meshgrid(1:height,1:width);
%     xx = vec(xx');
%     yy = vec(yy');
    quiver(xx2(:),yy2(:),uu(:),vv(:),4,'LineWidth',1, 'Color','r','MarkerSize',10);axis image
end

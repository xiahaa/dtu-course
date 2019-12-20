function [x1,x2] = automatic_matching(im1,im2)
    points1 = detectHarrisFeatures(im1);
    points2 = detectHarrisFeatures(im2);

    im3 = cat(2,im1,im2);

    points1 = points1.selectStrongest(1000);
    points2 = points2.selectStrongest(1000);

    imshow(im3); hold on;

    [features1,validPoints1] = extractFeatures(im1,points1);
    [features2,validPoints2] = extractFeatures(im2,points2);

    indexPairs = matchFeatures(features1,features2);

    matchedPoints1 = validPoints1(indexPairs(:,1));
    matchedPoints2 = validPoints2(indexPairs(:,2));

    x1 = matchedPoints1.Location;
    x2 = matchedPoints2.Location;

    plot(x1(:,1),x1(:,2),'ro');
    plot(x2(:,1)+size(im1,2),x2(:,2),'go');

    shift = size(im1,2);
    cmap = lines(5);
    k = 1;
    for i = 1:size(x1,1)
        ptdraw = [x1(i,1), x1(i,2);
                  x2(i,1)+shift, x2(i,2)];
        plot(ptdraw(:,1),ptdraw(:,2),'LineStyle','-','LineWidth',1,'Color',cmap(k,:));
        k = mod(k+1,5);if k == 0 k = 1;end
    end

    [fRANSAC, inliers] = estimateFundamentalMatrix(matchedPoints1,...
        matchedPoints2,'Method','RANSAC',...
        'NumTrials',5000,'DistanceThreshold',1e-2);

    inlierPoints1 = matchedPoints1(inliers);
    inlierPoints2 = matchedPoints2(inliers);

    figure
    imshow(im3); hold on;

    x1 = inlierPoints1.Location;
    x2 = inlierPoints2.Location;
    
end
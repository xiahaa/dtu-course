function test_on_knn_search
    %% setup vl_feat
    vl_setup;
    vl_version;
    
    %% sample
    X = rand(3,10000);
    kdtree = vl_kdtreebuild(X);
    
    Q = rand(3,1);
    tic
    [index, distance] = vl_kdtreequery(kdtree,X,Q,10);
    toc
    index
    distance
    tic
    dist = repmat(Q,1,size(X,2)) - X;
    dist = dist(1,:).^2+dist(2,:).^2+dist(3,:).^2;
    [minid,minval] = min(dist,10);
    toc
    minid
    minval
end
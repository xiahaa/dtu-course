function curve = suppressSelfIntersection(curve)
    cnt = size(curve,2);
%     [~,~,segments]=selfintersect(curve(2,:),curve(1,:)); 
    segments = checkSelfIntersect(curve(1,:),curve(2,:));
    id = ones(1,size(curve,2));
    for i = 1:size(segments,1)
        id(segments(i,1)+2:segments(i,2)-2) = 0;
    end
    id = id == 1;
    curve = curve(:,id);
    curve = reInterpolate(curve,cnt);
end

function segments = checkSelfIntersect(x,y)
    minx = min(x);
    maxx = max(x);
    miny = min(y);
    maxy = max(y);
    
    % bounding box
    bb = zeros(round(maxx - minx)+1,round(maxy-miny) + 1);
    
    % inter distance
    dist = sqrt((x(2:end) - x(1:end-1)).^2 + (y(2:end) - y(1:end-1)).^2);
    
    % 
    interpolateCnt = round(dist);
    
    % interpolation
    xf = zeros(sum(interpolateCnt),1);
    yf = zeros(sum(interpolateCnt),1);
    id = zeros(sum(interpolateCnt),1);
    k = 1;
    for i = 1:size(x,2)-1
        x1 = linspace(x(i),x(i+1),interpolateCnt(i));
        y1 = interp1([x(i);x(i+1)],[y(i);y(i+1)],x1);
        xf(k:k+length(x1)-1) = x1;
        yf(k:k+length(x1)-1) = y1;
        id(k:k+length(x1)-1) = i;
        k = k + length(x1);
    end
    
    %
    xf = round(xf-minx)+1; yf = round(yf-miny)+1;
    ind = round((yf-1).*size(bb,1)+xf);
    
    bb(ind) = id;
    idnew = bb(ind);
    
    intersects = abs(idnew - id) > 5 & abs(idnew - id) < 0.9*max(id);
        
    segments = [id(intersects) idnew(intersects)];
end
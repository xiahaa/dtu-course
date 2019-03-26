function curve = suppressSelfIntersection(curve)
    cnt = size(curve,2);
%     [~,~,segments]=selfintersect(curve(2,:),curve(1,:)); 
    segments = checkSelfIntersect(curve(1,:),curve(2,:));
    id = ones(1,size(curve,2));
    for i = 1:size(segments,1)
        id(segments(i,1):segments(i,2)) = 0;
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
    
    mask = poly2mask(y-miny+1,x-minx+1,size(bb,1),size(bb,2));
    boundary = findBoundary(mask);
    
    xx = round(x-minx)+1;
    yy = round(y-miny)+1;
    xx = xx';
    yy = yy';
    nnu = [xx-1, yy];n1 = nnu(:,1) > 0;
    nnd = [xx+1, yy];n2 = nnd(:,1) <= size(bb,1);
    nnl = [xx, yy-1];n3 = nnl(:,2) > 0;
    nnr = [xx, yy+1];n4 = nnr(:,2) <= size(bb,2);
    
    valid = zeros(size(x,2),1);valid = valid == 1;
    valid(n1) = valid(n1) | (boundary((nnu(n1,2)-1).*size(bb,1)+(nnu(n1,1))) == 1);
    valid(n2) = valid(n2) | (boundary((nnd(n2,2)-1).*size(bb,1)+(nnd(n2,1))) == 1);
    valid(n3) = valid(n3) | (boundary((nnl(n3,2)-1).*size(bb,1)+(nnl(n3,1))) == 1);
    valid(n4) = valid(n4) | (boundary((nnr(n4,2)-1).*size(bb,1)+(nnr(n4,1))) == 1);
    
    id = 1:size(x,2);
    idd = id(~valid);
    segments = [];
    if numel(idd) == 0
        return;
    end
    
    if numel(idd) > 0 && mod(numel(idd),2) == 0
        segments = [idd(1:2:end-1) idd(2:2:end)];
    else
        findidd = idd(end);
%         nextidd = min(findidd + 1,size(x,2));
        idd = [idd findidd];
        segments = [idd(1:2:end-1) idd(2:2:end)];
    end
        
    
end
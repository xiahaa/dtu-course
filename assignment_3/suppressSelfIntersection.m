function curve = suppressSelfIntersection(curve,varargin)
    if nargin == 3
        m = varargin{1};
        n = varargin{2};
    end

    cnt = size(curve,2);
%     [~,~,segments]=selfintersect(curve(2,:),curve(1,:)); 
    segments = InterX(curve,m,n);
    id = ones(1,size(curve,2));
    
    for i = 1:size(segments,1)
        id(segments(i,1):segments(i,2)) = 0;
    end
    id = id == 1;
    curve = curve(:,id);
    curve = reInterpolate(curve,cnt);
end

function segments = InterX(curve,m,n)
    point = [];
    id1 = [];
    for i = 1:size(curve,2)-1
        x = [curve(2,i) curve(2,i+1)];
        y = [curve(1,i) curve(1,i+1)];
        
        nPoints = max(abs(diff(x)), abs(diff(y)))+1;  % Number of points in line
        rIndex = round(linspace(y(1), y(2), nPoints));  % Row indices
        cIndex = round(linspace(x(1), x(2), nPoints));  % Column indices
        
        rcIndex = (cIndex-1).*m + rIndex;
        
        [rcIndex,id] = unique(rcIndex);
        rIndex = rIndex(id);
        cIndex = cIndex(id);
        
        if ~isempty(point) && rcIndex(1) == ((point(end,2)-1)*m+point(end,1))
            rIndex(1) = [];
            cIndex(1) = [];
        end
        point = [point;[rIndex' cIndex']];
        id1 = [id1;i*ones(length(rIndex),1)];
    end
    point(end,:) = [];
    id1(end) = [];
    
    bb = zeros(m,n);
    ind = round((point(:,2)-1).*m+point(:,1));
    bb(ind) = id1;
    idnew = bb(ind);
    intersects = abs(idnew - id1) > 5 & abs(idnew - id1) < 0.5*size(curve,2);
    
    segments = [id1(intersects) idnew(intersects)];
end
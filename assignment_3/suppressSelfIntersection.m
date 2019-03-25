function curve = suppressSelfIntersection(curve)
    [~,~,segments]=selfintersect(curve(2,:),curve(1,:)); 
    id = ones(1,size(curve,2));
    for i = 1:size(segments,1)
        id(segments(i,1)+1:segments(i,2)-1) = 0;
    end
    id = id == 1;
    curve = curve(:,id);
    curve = reInterpolate(curve);
end
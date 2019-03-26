function curve = reInterpolate(curve,Num)
    accumulateDist = zeros(1,size(curve,2));
    err = curve(:,[2:end]) - curve(:,1:end-1);
    dist = diag(err'*err)';
    for i = 2:size(curve,2)
        accumulateDist(i) = accumulateDist(i-1) + cumsum(dist(i-1));
    end
    accumulateDistNew = linspace(0,accumulateDist(end),Num);
    try
        xnew = interp1(accumulateDist,curve(1,:),accumulateDistNew);
    catch
%         error('s');
    end
    ynew = interp1(accumulateDist,curve(2,:),accumulateDistNew);
    curve = [xnew;ynew];
end
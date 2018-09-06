function currt = findVirtualTime(t, current_state, trajCoeff, trajTime, mode)

persistent sampleMode X ts total_time

if numel(X) == 0 | numel(ts) == 0 | numel(current_state) == 0
   total_time = trajTime(end);
   ts = trajTime;
   X = trajCoeff;
   sampleMode = mode;
   return
end

numCoeff = 6;
if sampleMode == 1
    currt = t;
    return;
elseif sampleMode == 2
    P = [ 91.1325   -0.0000   -0.0000    9.5331   -0.0000   -0.0000;
           -0.0000  128.5284   -0.0000   -0.0000   10.3061   -0.0000;
           -0.0000   -0.0000  111.8259   -0.0000   -0.0000    6.8330;
            9.5331   -0.0000   -0.0000   12.2852   -0.0000   -0.0000;
           -0.0000   10.3061   -0.0000   -0.0000   13.1432   -0.0000;
           -0.0000   -0.0000    6.8330   -0.0000   -0.0000    6.3460];

%     P = eye(6);
    
    persistent coarseSamplesP coarseSamplesT
    if numel(coarseSamplesP) == 0 || numel(coarseSamplesT) == 0
        tss = ts(1):0.01:ts(end);tss = tss';
        for j = 1:numel(tss)
            k = find(ts<=tss(j),1,'last');
            k(k==numel(ts)) = k(k==numel(ts)) - 1;
            scalar = 1./(ts(k+1)-ts(k));
            tnorm = (tss(j) - ts(k))*scalar;
            
            polyx = X((k-1)*numCoeff+1:k*numCoeff,1);
            polyy = X((k-1)*numCoeff+1:k*numCoeff,2);
            polyz = X((k-1)*numCoeff+1:k*numCoeff,3);
            
            xsamples = polyx(1) + polyx(2).*tnorm + polyx(3).*tnorm.^2 + polyx(4).*tnorm.^3 + polyx(5).*tnorm.^4 + polyx(6).*tnorm.^5;
            ysamples = polyy(1) + polyy(2).*tnorm + polyy(3).*tnorm.^2 + polyy(4).*tnorm.^3 + polyy(5).*tnorm.^4 + polyy(6).*tnorm.^5;
            zsamples = polyz(1) + polyz(2).*tnorm + polyz(3).*tnorm.^2 + polyz(4).*tnorm.^3 + polyz(5).*tnorm.^4 + polyz(6).*tnorm.^5;
            
            vxsamples = polyx(2) + polyx(3).*tnorm.*2 + polyx(4).*3.*tnorm.^2 + polyx(5).*4.*tnorm.^3 + polyx(6).*5.*tnorm.^4;
            vysamples = polyy(2) + polyy(3).*tnorm.*2 + polyy(4).*3.*tnorm.^2 + polyy(5).*4.*tnorm.^3 + polyy(6).*5.*tnorm.^4;
            vzsamples = polyz(2) + polyz(3).*tnorm.*2 + polyz(4).*3.*tnorm.^2 + polyz(5).*4.*tnorm.^3 + polyz(6).*5.*tnorm.^4;
            
            vxsamples = vxsamples.*scalar;
            vysamples = vysamples.*scalar;
            vzsamples = vzsamples.*scalar;
            
            coarseSamplesP(j,:) = [xsamples ysamples zsamples vxsamples vysamples vzsamples];
            coarseSamplesT(j) = tss(j);
        end
    end
    
    cs = current_state(1:6)';
    
    err = (repmat(cs, size(coarseSamplesP,1),1) - coarseSamplesP);
    dist = diag(err * P * err');
    [~,minid] = min(dist);
    
    tv = coarseSamplesT(minid);
    currt = tv;
elseif sampleMode == 3
    persistent sdist st sc
    if numel(sdist) == 0
        sc = 0;
        tss = ts(1):0.01:ts(end);
        oldp = [0 0 0];
        for j = 1:numel(tss)
            k = find(ts<=tss(j),1,'last');
            k(k==numel(ts)) = k(k==numel(ts)) - 1;
            scalar = 1./(ts(k+1)-ts(k));
            tnorm = (tss(j) - ts(k))*scalar;
            
            polyx = X((k-1)*numCoeff+1:k*numCoeff,1);
            polyy = X((k-1)*numCoeff+1:k*numCoeff,2);
            polyz = X((k-1)*numCoeff+1:k*numCoeff,3);
            
            xsamples = polyx(1) + polyx(2).*tnorm + polyx(3).*tnorm.^2 + polyx(4).*tnorm.^3 + polyx(5).*tnorm.^4 + polyx(6).*tnorm.^5;
            ysamples = polyy(1) + polyy(2).*tnorm + polyy(3).*tnorm.^2 + polyy(4).*tnorm.^3 + polyy(5).*tnorm.^4 + polyy(6).*tnorm.^5;
            zsamples = polyz(1) + polyz(2).*tnorm + polyz(3).*tnorm.^2 + polyz(4).*tnorm.^3 + polyz(5).*tnorm.^4 + polyz(6).*tnorm.^5;
            
            if j == 1
                sdist(j) = 0;
                st(j) = tss(j);
                oldp = [xsamples, ysamples, zsamples];
            else
                newp = [xsamples, ysamples, zsamples];
                sdist(j) = sdist(j-1) + norm(newp - oldp);
                st(j) = tss(j);
                oldp = newp;
            end
        end
    end
    
    cs = current_state(1:6)';
    preds = sc + norm(cs(4:6))*0.05;
    err = abs(sdist - preds);
    [~,minid] = min(err);
    tv = st(minid+5);
    currt = tv;
    sc = sdist(minid+5);
%     currt
end

end
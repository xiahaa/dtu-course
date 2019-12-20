function test_on_dop
    ii = 1;
    total = 1;
    success = zeros(total,1);
    den = zeros(total,1);
    for ii = 1:total
        n = round(10 + rand(1) * 30);
        m = round(4 + rand(1) * 6);
        
        p2 = rand(n, 3);
        %% normalize
        p2 = p2./sqrt(p2(:,1).^2+p2(:,2).^2+p2(:,3).^2);
        v1 = 1:1:n;
        C = nchoosek(v1,m);
        minGDOP = 1e6;
        minid = -1;
        maxGDOP2 = -1e6;
        maxid2 = -1;
%         tic
        for i = 1:size(C,1)
            A = p2(C(i,:)',:);
            GDOP = trace(inv(A'*A));
            if minGDOP > GDOP
                minGDOP = GDOP;
                minid = i;
            end
            if maxGDOP2 < det(A'*A)
                maxGDOP2 = det(A'*A);
                maxid2 = i;
            end
        end
%         toc
        o1 = C(minid,:);
        o2 = C(maxid2,:);

        %% idea
%         tic
        correlation = 1e6.*ones(n,n);
        for i = 1:n-1
            u1 = p2(i,:);
            urest = p2(i+1:end,:)';
            correlation(i,i+1:end) = u1*urest;
            correlation(i+1:end,i) = correlation(i,i+1:end)';
        end
        sum_correlation = sum(correlation,2) - diag(correlation);
        [val,id] = min(sum_correlation);
        sum_correlation(id) = 1e8;
        indices = [id];
        cnt = 1;
        while cnt < m
            ids = 1:n;
%             ids(id) = [];
            sum_correlation = sum_correlation;% - correlation(ids,id);
            [val,id] = min(sum_correlation);
            sum_correlation(id) = 1e8;
            indices = [indices;id];
            cnt = cnt + 1;
        end
%         o2 = indices';
        o1 = sort(o1);
        o2 = sort(o2);
        success(ii) = numel(find(abs(o2-o1)<1));
        den(ii) = numel(o1);
%         toc
    end
    sum(success) / sum(den)
end
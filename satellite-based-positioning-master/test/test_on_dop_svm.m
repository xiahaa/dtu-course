function test_on_dop_svm
    n = 30;
    dim = 3;
    m = 5;
    runcnt = 10000;
    k = 1;
    rs = zeros(runcnt, n);
    
    bins1 = 72; bins2 = 72; bins3 = 4;
    
    fea_vecs1 = zeros(runcnt, n, bins1);
    fea_vecs2 = zeros(runcnt, n, bins2);
    fea_vecs3 = zeros(runcnt, n, bins3);
    
    parfor k = 1:runcnt
        random_vecs = gen_ran_vecs(n, dim);
        p = v_normalize(random_vecs);
        r = rank_sat(p, m);
        rs(k,:) = r';
        fea_vec1 = extrac_feature_vecs1(p, bins1);
        fea_vecs1(k,:,:) = fea_vec1;
        fea_vec2 = extrac_feature_vecs2(p, bins2);
        fea_vecs2(k,:,:) = fea_vec2;
        fea_vec3 = extrac_feature_vecs3(p, bins3);
        fea_vecs3(k,:,:) = fea_vec3;
    end
    
    save('trainingdata.mat','rs','fea_vecs1','fea_vecs2','fea_vecs3');
    
    %% save
    currpath = pwd();
    
    fid = fopen(strcat(currpath,'/training_f1.dat'),'w');
    for k = 1:runcnt
        r = rs(k,:)';
        fea_vec = fea_vecs1(k,:,:);
        fea_vec = reshape(fea_vec,n,bins1,1);
        for i = 1:numel(r)
            fprintf(fid,'%d qid:%d ', r(i), k);
            for j = 1:size(fea_vec,2)
                fprintf(fid,'%d:%f ', j, fea_vec(i,j));
            end
            fprintf(fid,'# %d\n', i);
        end
    end
    fclose(fid);
    
    fid = fopen(strcat(currpath,'/training_f2.dat'),'w');
    for k = 1:runcnt
        r = rs(k,:)';
        fea_vec = fea_vecs2(k,:,:);
        fea_vec = reshape(fea_vec,n,bins2,1);
        for i = 1:numel(r)
            fprintf(fid,'%d qid:%d ', r(i), k);
            for j = 1:size(fea_vec,2)
                fprintf(fid,'%d:%f ', j, fea_vec(i,j));
            end
            fprintf(fid,'# %d\n', i);
        end
    end
    fclose(fid);
    
    fid = fopen(strcat(currpath,'/training_f3.dat'),'w');
    for k = 1:runcnt
        r = rs(k,:)';
        fea_vec = fea_vecs3(k,:,:);
        fea_vec = reshape(fea_vec,n,bins3,1);
        for i = 1:numel(r)
            fprintf(fid,'%d qid:%d ', r(i), k);
            for j = 1:size(fea_vec,2)
                fprintf(fid,'%d:%f ', j, fea_vec(i,j));
            end
            fprintf(fid,'# %d\n', i);
        end
    end
    fclose(fid);
    
end

function random_vecs = gen_ran_vecs(n, dim)
    rng('shuffle');
    random_vecs = rand(n,dim) * 2 - 1;
end

function nvecs = v_normalize(vecs)
    nvecs = vecs./vecnorm(vecs,2,2);
end

function r = rank_sat(p, m)
    n = size(p,1);
    v1 = 1:1:n;
    C = nchoosek(v1,m);
    
    %% find optimum
    minGDOP = 1e6;
    minid = -1;
    for i = 1:size(C,1)
        A = p(C(i,:)',:);
        GDOP = trace(inv(A'*A));
        if minGDOP > GDOP
            minGDOP = GDOP;
            minid = i;
        end
    end
    opt = C(minid,:);
    
    r = zeros(n,1);
    r(opt,:) = n - m + 1;
    
    %% incrementally compute GDOP
    indices = opt;
    restindices = 1:n;
    restindices(indices) = [];
    for i = 1:(n-m)
        newindices = [repmat(indices, numel(restindices), 1) restindices'];
        minGDOP = 1e6;
        minid = -1;
        for j = 1:numel(restindices)
            A = p(newindices(j,:),:);
            GDOP = trace(inv(A'*A));
            if minGDOP > GDOP
                minGDOP = GDOP;
                minid = j;
            end
        end
        indices = [indices restindices(minid)];
        r(restindices(minid)) = n - m + 1 - i;
        restindices(minid) = [];
    end
end

function fea_vecs = extrac_feature_vecs2(p, bins)
    fea_vecs = zeros(size(p,1),bins);
    bin = linspace(-1,1,bins + 1);
    bin = bin(2:end);
    n = size(p,1);
    for i = 1:n
        u1 = p(i,:);
        ids = 1:n;
        ids(i) = [];
        urest = p(ids,:)';
        correlation = u1*urest;
        correlation = sort(correlation,'ascend');
        fea_vecs(i,1:numel(correlation)) = correlation;
    end
    %% normalization
%     fea_vecs = fea_vecs ./ (n-1);
end

function fea_vecs = extrac_feature_vecs3(p, bins)
    fea_vecs = zeros(size(p,1),bins);
    bin = linspace(-1,1,bins + 1);
    bin = bin(2:end);
    n = size(p,1);
    for i = 1:n
        u1 = p(i,:);
        ids = 1:n;
        ids(i) = [];
        urest = p(ids,:)';
        correlation = u1*urest;
        
        f1 = sum(correlation);
        f2 = std(correlation);
        f3 = max(correlation);
        f4 = min(correlation);
        
        fea_vecs(i,1:4) = [f1 f2 f3 f4];
    end
    %% normalization
%     fea_vecs = fea_vecs ./ (n-1);
end

function fea_vecs = extrac_feature_vecs1(p, bins)
    fea_vecs = zeros(size(p,1),bins);
    bin = linspace(-1,1,bins + 1);
    bin = bin(2:end);
    n = size(p,1);
    for i = 1:n
        u1 = p(i,:);
        ids = 1:n;
        ids(i) = [];
        urest = p(ids,:)';
        correlation = u1*urest;
        
        for j = 1:numel(correlation)
            bin_id = find(bin >= correlation(j),1);
            fea_vecs(i,bin_id) = fea_vecs(i,bin_id) + 1;
        end
    end
    %% normalization
    fea_vecs = fea_vecs ./ (n-1);
end
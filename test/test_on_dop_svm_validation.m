function test_on_dop_svm_validation
    n = 20;
    dim = 3;
    m = 8;
    bins = 4;
    runcnt = 100;
    k = 1;
    
    bins1 = 72; bins2 = 72; bins3 = 4;
    
    clsy1 = zeros(bins1,1);
    clsy2 = zeros(bins2,1);
    clsy3 = zeros(bins3,1);
        
    succf1 = zeros(runcnt, 1);
    succf2 = zeros(runcnt, 1);
    succf3 = zeros(runcnt, 1);
    
    succ3 = zeros(runcnt, 1);
    succ4 = zeros(runcnt, 1);
    
    total = zeros(runcnt, 1);
        
    clsy1 = read_model('./svm_rank_windows/model_f1.dat',bins1); 
    clsy2 = read_model('./svm_rank_windows/model_f2.dat',bins2); 
    clsy3 = read_model('./svm_rank_windows/model_f3.dat',bins3); 
    
    for k = 1:runcnt
        random_vecs = gen_ran_vecs(n, dim);
        p = v_normalize(random_vecs);
        fea_vec1 = extrac_feature_vecs(p, bins1);
        fea_vec2 = extrac_feature_vecs2(p, bins2);
        fea_vec3 = extrac_feature_vecs3(p, bins3);
        [sf1,sf2,sf3,s2,s3,t] = validation(p, m, fea_vec1, clsy1, fea_vec2, clsy2, fea_vec3, clsy3);
        
        succf1(k) = sf1;succf2(k) = sf2;succf3(k) = sf3;
        succ3(k) = s2;
        succ4(k) = s3;
        total(k) = t;
    end
    
   sum(succf1) / sum(total)
   sum(succf2) / sum(total)
   sum(succf3) / sum(total)
   sum(succ3) / sum(total)
   sum(succ4) / sum(total)
end

function clsy = read_model(file, bins)
    fid = fopen(file,'r');
    row = 12;
    currow = 1;
    clsy = zeros(bins,1);
    while ~feof(fid)
        tline = fgetl(fid);
        if currow ~= row
            currow = currow + 1;
            continue;
        end
        j = 3;
        while j <= size(tline,2)-2
            line = '';
            while tline(1,j) ~= ' '
                line = strcat(line,tline(1,j));
                j = j + 1;
            end
            res = strsplit(line,':');
            id = str2num(res{1});
            w = str2double(res{2});
            clsy(id) = w;
            j = j + 1;
        end
    end
end

function random_vecs = gen_ran_vecs(n, dim)
    rng('shuffle');
    random_vecs = rand(n,dim) * 2 - 1;
end

function nvecs = v_normalize(vecs)
    nvecs = vecs./vecnorm(vecs,2,2);
end

function [sf1,sf2,sf3,s2,s3,t] = validation(p, m, fea_vec1, clsy1, fea_vec2, clsy2, fea_vec3, clsy3)
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

    scores = fea_vec1 * clsy1;
    [val, id1] = sort(scores,'descend');
    opt1 = id1(1:m);
    
    scores = fea_vec2 * clsy2;
    [val, id2] = sort(scores,'descend');
    opt2 = id2(1:m);
    
    scores = fea_vec3 * clsy3;
    [val, id3] = sort(scores,'descend');
    opt3 = id3(1:m);
    
    correlation1 = 1e6.*ones(n,n);
    for i = 1:n-1
        u1 = p(i,:);
        urest = p(i+1:end,:)';
        correlation1(i,i+1:end) = (u1*urest);
        correlation1(i+1:end,i) = correlation1(i,i+1:end)';
    end
    sum_correlation1 = sum(correlation1,2) - diag(correlation1);
    [val,id4] = sort(sum_correlation1,'ascend');
    opt4 = id4(1:m)';
    
    correlation = 1e6.*ones(n,n);
    for i = 1:n-1
        u1 = p(i,:);
        urest = p(i+1:end,:)';
        correlation(i,i+1:end) = 2.*(u1*urest)-1;
        correlation(i+1:end,i) = correlation(i,i+1:end)';
    end
    sum_correlation = sum(correlation,2) - diag(correlation);
    [val,id5] = sort(sum_correlation,'ascend');
    opt5 = id5(1:m)';
    
    sf1 = 0;sf2 = 0;sf3 = 0;
    s2 = 0;
    s3 = 0;
    for i = 1:m
        if ~isempty(find(abs(opt-opt1(i))<1))
            sf1 = sf1 + 1;
        end
        if ~isempty(find(abs(opt-opt2(i))<1))
            sf2 = sf2 + 1;
        end
        if ~isempty(find(abs(opt-opt3(i))<1))
            sf3 = sf3 + 1;
        end
        if ~isempty(find(abs(opt-opt4(i))<1))
            s2 = s2 + 1;
        end
        if ~isempty(find(abs(opt-opt5(i))<1))
            s3 = s3 + 1;
        end
    end
    
%     s1 = numel(find(abs(opt2-opt)<1));
%     s2 = numel(find(abs(opt3-opt)<1));
    t = numel(opt);
    
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

function fea_vecs = extrac_feature_vecs(p, bins)
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




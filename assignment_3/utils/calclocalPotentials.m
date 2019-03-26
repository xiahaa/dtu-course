function localPotential = calclocalPotentials(im, configuration, f, miuf, alpha, beta)
%Compute potentials per pixels for each label.
    localPotential = zeros(size(im,1),size(im,2),numel(f));
    % likelihood
    for i = 1:numel(f)
        p1(:,:,i) = alpha.*(miuf(f(i))-double(im)).^2;
    end
    [m,n] = size(im);
    % smoothness
    [xx,yy]=meshgrid(1:m,1:n);
    xx = vec(xx'); yy = vec(yy');
    nnu = [xx-1, yy]; n1 = nnu(:,1) > 0;
    nnd = [xx+1, yy]; n2 = nnd(:,1) <= m;
    nnl = [xx, yy-1]; n3 = nnl(:,2) > 0;
    nnr = [xx, yy+1]; n4 = nnr(:,2) <= n;
    
    localPotential = p1;
    
    if ~isempty(configuration)
        for i = 1:numel(f)
            b = zeros(m*n,1);
        %     fi = configuration((yy-1).*m+xx);
            b(n1) = b(n1) + double(configuration((nnu(n1,2)-1).*m + nnu(n1,1)) ~= f(i));
            b(n2) = b(n2) + double(configuration((nnd(n2,2)-1).*m + nnd(n2,1)) ~= f(i));
            b(n3) = b(n3) + double(configuration((nnl(n3,2)-1).*m + nnl(n3,1)) ~= f(i));
            b(n4) = b(n4) + double(configuration((nnr(n4,2)-1).*m + nnr(n4,1)) ~= f(i));
            b = reshape(b,m,n);

            localPotential(:,:,i) = localPotential(:,:,i) + b.*beta;
        end
    end
end
function boundary = findBoundary(bw)
%finding boundary.
    boundary = zeros(size(bw));
    m = size(bw,1);
    n = size(bw,2);
    % mask pixels
    [row,col] = find(bw);
    nnu = [row-1,col];n1 = nnu(:,1) > 0;
    nnd = [row+1,col];n2 = nnd(:,1) <= m;
    nnl = [row,col-1];n3 = nnl(:,2) > 0;
    nnr = [row,col+1];n4 = nnr(:,2) <= n;
    % 4 neighbors
    valid = zeros(numel(row),1);valid = valid == 1;
    valid(n1) = valid(n1) | bw((nnu(n1,2)-1).*m+nnu(n1,1)) == 0;
    valid(n2) = valid(n2) | bw((nnd(n2,2)-1).*m+nnd(n2,1)) == 0;
    valid(n3) = valid(n3) | bw((nnl(n3,2)-1).*m+nnl(n3,1)) == 0;
    valid(n4) = valid(n4) | bw((nnr(n4,2)-1).*m+nnr(n4,1)) == 0;
    % to logical
    boundary(col(valid).*m+row(valid)) = 1;
end
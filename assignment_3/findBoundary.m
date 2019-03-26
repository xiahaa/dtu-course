% poisson blending
function boundary = findBoundary(bw)
    boundary = zeros(size(bw));
    for i = 2:size(bw,1)-1
        for j = 2:size(bw,2)-1
            block = bw(i-1:i+1,j-1:j+1);
            if block(5) == 1 && sum(block([1,2,3,4,6,7,8,9])==0)~=0
                boundary(i,j) = 1;
            end
        end
    end
end
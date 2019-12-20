function TV = total_variantion_calculator(I)
% compute the total variantion of an image using vectorization.
% Author: xiahaa@space.dtu.dk
    [M,N] = size(I);
    D1 = spdiags([-ones(M-1,1),ones(M-1,1);0 1],[0,1],M,M);
    D2 = spdiags([-ones(N-1,1),ones(N-1,1);0 1],[0,1],N,N);
    D = vertcat(kron(speye(N),D1),kron(D2,speye(M)));
    ftv = @(x) sum(sqrt(sum(reshape(D*x,[],2).^2)));
    TV = ftv(I(:));
end
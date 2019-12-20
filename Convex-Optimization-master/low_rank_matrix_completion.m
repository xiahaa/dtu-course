function low_rank_matrix_completion
    rng(234923);    % for reproducible results
    N   = 16;       % the matrix is N x N
    r   = 2;        % the rank of the matrix
    df  = 2*N*r - r^2;  % degrees of freedom of a N x N rank r matrix
    nSamples    = 3*df; % number of observed entries

    % For this demo, we will use a matrix with integer entries
    % because it will make displaying the matrix easier.
    iMax    = 5;
    X       = randi(iMax,N,r)*randi(iMax,r,N); % Our target matrix
    
    rPerm   = randperm(N^2); % use "randsample" if you have the stats toolbox
    omega   = sort( rPerm(1:nSamples) );
    
    Y = nan(N);
    Y(omega) = X(omega);
    disp('The "NaN" entries represent unobserved values');
    disp(Y)
    
    observations = X(omega);    % the observed entries
    mu = .001;        % smoothing parameter

    % The solver runs in seconds
    tic
    Xk = solver_sNuclearBP( {N,N,omega}, observations, mu );
    toc
    disp('Recovered matrix (rounding to nearest .0001):')
    disp( round(Xk*10000)/10000 )
    % and for reference, here is the original matrix:
    disp('Original matrix:')
    disp( X )

    % The relative error (without the rounding) is quite low:
    fprintf('Relative error, no rounding: %.8f%%\n', norm(X-Xk,'fro')/norm(X,'fro')*100 );
% r = 2; % the rank;
% 
% N = 32; % the dimension
% 
% M = 32;
% 
% a = randn(N,r);
% 
% b = randn(M,r);
% 
% X = a*b'; % low rank Matrix;
% 
% A = rand(20,N);
% 
% Y = A*X;
% 
% % low rank approximation using nuclear norm
% 
% cvx_begin
% variable Xe(N,M)
% minimize norm_nuc(Xe)
% subject to 
%     A*Xe == Y;
% cvx_end
end
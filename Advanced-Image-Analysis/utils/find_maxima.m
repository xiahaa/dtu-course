function [varargout] = find_maxima(I, half_win_size)
%     I = abs(I);
    val = I(:);
    meanval = mean(val);
    stdv = var(val);
    thresh = meanval+2*stdv;%0.3*max(I(:));
    mx = imdilate(I, strel('square',2*half_win_size+1));

    % Make mask to exclude points within radius of the image boundary. 
    bordermask = zeros(size(I));
    bordermask(half_win_size+1:end-half_win_size, half_win_size+1:end-half_win_size) = 1;
    
    % Find maxima, threshold, and apply bordermask
    cimmx = (I==mx) & (I>thresh) & bordermask;
    [r,c] = find(cimmx);        % Find row,col coords.
    blobs_center = [r c];
    
    % visulization
    cc = blobs_center;
    num_segments = size(cc, 1);
    [M,N] = size(I);
    Idebug = cat(3,I,I,I);
    for i = 1:num_segments
        indices = sub2ind(size(I), cc(i,1), cc(i,2));
        nn = [-M +M -1 +1];% draw also 4 neighbors
        indices = [indices];% indices+nn(1) indices+nn(2) indices+nn(3) indices+nn(4)];
        indices = indices(indices > 0 & indices < M*N);
        Idebug(indices) = 1;
        Idebug(indices+M*N) = 0;
        Idebug(indices+M*N*2) = 0;
    end
    varargout{1} = blobs_center;
    if nargout == 2
        varargout{2} = Idebug;
    end
end
function [blobs_center] = detect_blobs(I, half_win_size)
    I = abs(I);
    val = I(:);
    meanval = mean(val);
    stdv = var(val);
    thresh = meanval+stdv;
    mx = imdilate(I, strel('square',2*half_win_size+1));

    % Make mask to exclude points within radius of the image boundary. 
    bordermask = zeros(size(I));
    bordermask(half_win_size+1:end-half_win_size, half_win_size+1:end-half_win_size) = 1;
    
    % Find maxima, threshold, and apply bordermask
    cimmx = (I==mx) & (I>thresh) & bordermask;
    [r,c] = find(cimmx);        % Find row,col coords.
    blobs_center = [r c];
end
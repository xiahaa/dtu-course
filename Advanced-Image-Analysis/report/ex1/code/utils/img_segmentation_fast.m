function varargout = img_segmentation_fast(I)
% run breadth-first search for contour finding and boundary length
% computation.
% Author: xiahaa@space.dtu.dk
    [M,N] = size(I);
    visited = zeros(M*N,1);
    current_label = 0;
    p = I(:);
    % 4 neightbor
    nnid = [-M,M,-1,+1];% left, right, up, down
    LIFO = zeros(M*N,1);
    LIFO_SIZE = 0;
    num_segments = 0;
    segments = {};
        
    for i = 1:M*N
        if visited(i) ~= 0 || p(i) == 0 continue; end
        current_label = current_label + 1;
        LIFO_SIZE = LIFO_SIZE + 1;LIFO(LIFO_SIZE)=i;
        buf = zeros(M*N,1);
        len = 0;
        
        while LIFO_SIZE>0
            ii = LIFO(1:LIFO_SIZE);LIFO_SIZE = 0;
            buf(len+1:len+numel(ii))=ii;len = len + numel(ii);
            visited(ii) = 1;
            nnl = ii + nnid(1);nnr = ii + nnid(2);nnu = ii + nnid(3);nnd = ii + nnid(4);
%             nn = sub2ind([M,N],ii) + nnid;
%             valid = nnl(nnl > 0 & nnl < M*N);
%             nn = nn(valid);
            nn = [nnl(nnl>0&nnl<M*N); ...
                  nnr(nnr>0&nnr<M*N); ...
                  nnu(nnu>0&nnu<M*N); ...
                  nnd(nnd>0&nnd<M*N)];
            nn = unique(nn')';
%             ii = LIFO(LIFO_SIZE);LIFO_SIZE = LIFO_SIZE-1;
%             len = len + 1; buf(len)=ii; visited(ii) = 1;
%             nn = ii + nnid;
%             nn = nn((nn > 0 & nn < M*N));
            nnvalid = nn(visited(nn) == 0 & p(nn) ~= 0);
            LIFO(LIFO_SIZE+1:LIFO_SIZE+numel(nnvalid))=nnvalid;
            LIFO_SIZE = LIFO_SIZE + numel(nnvalid);
        end
        % find one segment
        num_segments = num_segments + 1;
        [II,JJ] = ind2sub([M,N],buf(1:len));
        point_set = [II JJ];
        segment.point_set = point_set';
        segments{num_segments} = segment;
    end
%     
%     segments = {};
%     visited = zeros(size(I));
%     current_label = 0;
%     [M,N] = size(I);
%     LIFO = zeros(M*N,2);
%     LIFO_SIZE = 0;
%     num_segments = 0;
%     for i = 1:M
%         for j = 1:N
%             if visited(i,j) ~= 0 || I(i,j) == 0
%                 continue;
%             end
%             % not yet visited
%             n = 0;
%             segment.label = current_label;
%             point_set = [];% make this 2xN
%             segment.boundary_length = 0;
%             current_label = current_label + 1;
%             LIFO_SIZE=LIFO_SIZE+1;LIFO(LIFO_SIZE,:) = [i,j];
%             while LIFO_SIZE > 0
%                 iijj = LIFO(LIFO_SIZE,:);LIFO_SIZE = LIFO_SIZE - 1;
%                 ci = iijj(1);cj = iijj(2);
%                 % store
%                 n = n+1;
%                 point_set(n,:) = iijj;
%                 %
%                 visited(iijj(1),iijj(2)) = 1;
%                 % check neighboring
%                 nni = iijj(1)-1;
%                 nnj = iijj(2);
%                 % up
%                 if ( (nni > 0) && visited(nni, nnj) == 0 && I(nni,nnj) ~= 0)
%                     LIFO_SIZE=LIFO_SIZE+1;LIFO(LIFO_SIZE,:) = [nni,nnj];
%                 end
%                 % down
%                 nni = iijj(1)+1;
%                 nnj = iijj(2);
%                 if ( (nni <= M) && visited(nni, nnj) == 0 && I(nni,nnj) ~= 0)
%                     LIFO_SIZE=LIFO_SIZE+1;LIFO(LIFO_SIZE,:) = [nni,nnj];
%                 end
%                 % left
%                 nni = iijj(1);
%                 nnj = iijj(2)-1;
%                 if ( (nnj > 0) && visited(nni, nnj) == 0 && I(nni,nnj) ~= 0)
%                     LIFO_SIZE=LIFO_SIZE+1;LIFO(LIFO_SIZE,:) = [nni,nnj];
%                 end
%                 % right
%                 nni = iijj(1);
%                 nnj = iijj(2)+1;
%                 if ( (nnj <= N) && visited(nni, nnj) == 0 && I(nni,nnj) ~= 0)
%                     LIFO_SIZE=LIFO_SIZE+1;LIFO(LIFO_SIZE,:) = [nni,nnj];
%                 end
%             end
%             % find one segment
%             num_segments = num_segments + 1;
%             segment.point_set = point_set';
%             segments{num_segments} = segment;
%         end
%     end
    varargout{1} = segments;
end
function varargout = img_segmentation(I)
% run breadth-first search for contour finding and boundary length
% computation.
% Author: xiahaa@space.dtu.dk
    segments = {};
    visited = zeros(size(I));
    current_label = 0;
    [M,N] = size(I);
    LIFO = zeros(M*N,2);
    LIFO_SIZE = 0;
    num_segments = 0;
    for i = 1:M
        for j = 1:N
            if visited(i,j) ~= 0
                continue;
            end
            % not yet visited
            n = 0;
            segment.label = current_label;
            point_set = [];% make this 2xN
            segment.boundary_length = 0;
            current_label = current_label + 1;
            LIFO_SIZE=LIFO_SIZE+1;LIFO(LIFO_SIZE,:) = [i,j];
            while LIFO_SIZE > 0
                iijj = LIFO(LIFO_SIZE,:);LIFO_SIZE = LIFO_SIZE - 1;
                ci = iijj(1);cj = iijj(2);
                % store
                n = n+1;
                point_set(n,:) = iijj;
                %
                visited(iijj(1),iijj(2)) = 1;
                % check neighboring
                nni = iijj(1)-1;
                nnj = iijj(2);
                isbounder = 0;
                % up
                if ( (nni > 0) && visited(nni, nnj) == 0 )
                    if (I(nni, nnj) == I(ci,cj))
                        LIFO_SIZE=LIFO_SIZE+1;LIFO(LIFO_SIZE,:) = [nni,nnj];
                    else
                        isbounder = 1;
                    end
                end
                % down
                nni = iijj(1)+1;
                nnj = iijj(2);
                if ( (nni <= M) && visited(nni, nnj) == 0 )
                    if (I(nni, nnj) == I(ci,cj))
                        LIFO_SIZE=LIFO_SIZE+1;LIFO(LIFO_SIZE,:) = [nni,nnj];
                    else
                        isbounder = 1;
                    end
                end
                % left
                nni = iijj(1);
                nnj = iijj(2)-1;
                if ( (nnj > 0) && visited(nni, nnj) == 0 )
                    if (I(nni, nnj) == I(ci,cj))
                        LIFO_SIZE=LIFO_SIZE+1;LIFO(LIFO_SIZE,:) = [nni,nnj];
                    else
                        isbounder = 1;
                    end
                end
                % right
                nni = iijj(1);
                nnj = iijj(2)+1;
                if ( (nnj <= N) && visited(nni, nnj) == 0 )
                    if (I(nni, nnj) == I(ci,cj))
                        LIFO_SIZE=LIFO_SIZE+1;LIFO(LIFO_SIZE,:) = [nni,nnj];
                    else
                        isbounder = 1;
                    end
                end
                
                if isbounder == 1
                    segment.boundary_length = segment.boundary_length + 1;
                end
            end
            % find one segment
            num_segments = num_segments + 1;
            segment.point_set = point_set';
            segments{num_segments} = segment;
        end
    end
    varargout{1} = segments;
end
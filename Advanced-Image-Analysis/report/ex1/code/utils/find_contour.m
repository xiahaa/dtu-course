
function varargout = find_contour(Ib)
% run breadth-first search for contour finding.
% computation.
% Author: xiahaa@space.dtu.dk
    contours = {};
    visited = zeros(size(Ib));
    current_label = 0;
    [M,N] = size(Ib);
    num_contours = 0;
    LIFO = zeros(M*N,2);
    LIFO_SIZE = 0;
    for i = 1:M
        for j = 1:N
            if visited(i,j) ~= 0 || Ib(i,j) == 0
                continue;
            end
            % not yet visited
            n = 0;
            contour.label = current_label;
            point_set = [];% make this 2xN
            current_label = current_label + 1;
            LIFO_SIZE = LIFO_SIZE + 1;LIFO(LIFO_SIZE, :) = [i,j];
            
            while LIFO_SIZE ~= 0
                iijj = LIFO(LIFO_SIZE,:);LIFO_SIZE = LIFO_SIZE - 1;
                ci = iijj(1);cj = iijj(2);
                %
                visited(iijj(1),iijj(2)) = 1;
                % check neighboring
                nni = iijj(1)-1;
                nnj = iijj(2);
                isbounder = 0;
                % up
                if ( (nni > 0) && visited(nni, nnj) == 0 )
                    if (Ib(nni, nnj) == Ib(ci,cj))
                        LIFO_SIZE = LIFO_SIZE + 1;LIFO(LIFO_SIZE, :) = [nni,nnj];
                    else
                        isbounder = 1;
                    end
                end
                % down
                nni = iijj(1)+1;
                nnj = iijj(2);
                if ( (nni <= M) && visited(nni, nnj) == 0 )
                    if (Ib(nni, nnj) == Ib(ci,cj))
                        LIFO_SIZE = LIFO_SIZE + 1;LIFO(LIFO_SIZE, :) = [nni,nnj];
                    else
                        isbounder = 1;
                    end
                end
                % left
                nni = iijj(1);
                nnj = iijj(2)-1;
                if ( (nnj > 0) && visited(nni, nnj) == 0 )
                    if (Ib(nni, nnj) == Ib(ci,cj))
                        LIFO_SIZE = LIFO_SIZE + 1;LIFO(LIFO_SIZE, :) = [nni,nnj];
                    else
                        isbounder = 1;
                    end
                end
                % right
                nni = iijj(1);
                nnj = iijj(2)+1;
                if ( (nnj <= N) && visited(nni, nnj) == 0 )
                    if (Ib(nni, nnj) == Ib(ci,cj))
                        LIFO_SIZE = LIFO_SIZE + 1;LIFO(LIFO_SIZE, :) = [nni,nnj];
                    else
                        isbounder = 1;
                    end
                end
                
                if isbounder == 1
                    % store
                    n = n+1;
                    point_set(n,:) = iijj;
                end
            end
            
            if size(point_set,1) > 800 % minimum size for one solid contour
                % find one segment
                num_contours = num_contours + 1;
                contour.point_set = point_set';
                contours{num_contours} = contour;
            end
        end
    end
    varargout{1} = contours;
end

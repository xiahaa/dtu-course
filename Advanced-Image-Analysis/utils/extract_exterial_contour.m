
function varargout = extract_exterial_contour(contour)
    n = size(contour,2);
    contour = int32(contour);
    [~,sortid] = sort(contour(1,:),2);% sort by row
    contour = contour(:,sortid);
    exterial_contour = [];
    % first and last
    crow = contour(1,1);
    drow = contour(1,:) - crow;
    id = find(drow == 0);% same row
    exterial_contour = [exterial_contour contour(:,id)];
    i = 1 + numel(id);
    while i < n
        crow = contour(1,i);
        drow = contour(1,:) - crow;
        id = find(drow == 0);% same row
        if numel(id) == 1
            exterial_contour = [exterial_contour contour(:,[id])];
            continue;
        end
        [~,minid] = min(contour(2,id));
        [~,maxid] = max(contour(2,id));
        exterial_contour = [exterial_contour contour(:,[id(minid),id(maxid)])];
        i = i + numel(id);
    end
    crow = contour(1,end);
    drow = contour(1,:) - crow;
    id = find(drow == 0);% same row
    exterial_contour = [exterial_contour contour(:,id)];
    varargout{1} = double(exterial_contour);
end
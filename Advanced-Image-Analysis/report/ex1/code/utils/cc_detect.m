
function varargout = cc_detect(I)
    % raw image to binary image
    Ib = imbinarize(I);
    % find contour
    contours = find_contour(Ib);
    % visulization, only for debuging
    num_contours = size(contours, 2);
    valid_contours = [];
    ss = [];
    cc = [];
    for i = 1:num_contours
        contour = contours{i};
        contour.point_set = extract_exterial_contour(contour.point_set);
%         scatter(contour.point_set(2,:),contour.point_set(1,:));hold on;
        [e_params, ~] = fit_ellipse(contour.point_set);
        [canonical_params,~] = AtoG(e_params);
        semi_a = canonical_params(3);
        semi_b = canonical_params(4);
        if abs(semi_b/semi_a - 1) < 0.2
            valid_contours = [valid_contours;i];
            ss = [ss;abs(semi_b/semi_a - 1)];
            cc = [cc;[canonical_params']];
        end
%         scatter(canonical_params(2,:),canonical_params(1,:));
%         plotellipse(canonical_params(1),canonical_params(2),canonical_params(3),canonical_params(4),canonical_params(5));
    end
    num_contours = size(valid_contours,1);
    if num_contours > 1
        [~,minid] = min(ss);
    else
        minid = 1;
    end
    valid_contours = valid_contours(minid);
    num_contours = size(valid_contours,1);
    % imshow(Ic,[]);
%     ccmap = jet(num_contours);
%     Ic = cat(3, Ib, Ib, Ib);
%     Ic = im2double(Ic);
%     [M,N] = size(Ib);
%     for i = 1:num_contours
%         contour = contours{valid_contours(i)};
%         indices = sub2ind(size(Ib), contour.point_set(1,:), contour.point_set(2,:));
%         Ic(indices) = ccmap(i,1);
%         Ic(indices+M*N) = ccmap(i,2);
%         Ic(indices+M*N*2) = ccmap(i,3);
%     end
%     imshow(Ic,[]);hold on;
%     plot(cc(minid,2),cc(minid,1),'rx','MarkerSize',10);
    varargout{1} = cc(minid,:);
end
function normals = computeNormal(curve)
%Naive method of computing normals
	% firstly, computing tangent direction
    tangent = zeros(2,size(curve,2));
    tangent(:,1) = curve(:,2) - curve(:,end);
    tangent(:,2:end-1) = curve(:,3:end) - curve(:,1:end-2);
    tangent(:,end) = curve(:,1) - curve(:,end-1);
    % normalize
    tangent = tangent ./ sqrt(tangent(1,:).^2+tangent(2,:).^2);
    % rotate to outward normal direction
    normals = [tangent(2,:);-tangent(1,:)];
end
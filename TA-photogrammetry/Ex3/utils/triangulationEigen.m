function P = triangulationEigen(x1,P1,x2,P2,varargin)
%Implementation of Linear-Eigen triangulation method.
    if nargin > 4
        w1 = varargin{1};
        w2 = varargin{2};
    else
        w1 = 1;
        w2 = 1;
    end

    A(1,:) = (x1(1).*P1(3,:) - P1(1,:)).*w1;
    A(2,:) = (x1(2).*P1(3,:) - P1(2,:)).*w1;
    A(3,:) = (x2(1).*P2(3,:) - P2(1,:)).*w2;
    A(4,:) = (x2(2).*P2(3,:) - P2(2,:)).*w2;
    [~,~,V] = svd(A);
    P = V(:,end);
    P = P(1:3)./P(end);
end
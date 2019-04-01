
function P = triangulationLS(x1,P1,x2,P2,varargin)
%Implementation of Linear-least-square triangulation method.
    if nargin > 4
        w1 = varargin{1};
        w2 = varargin{2};
    else
        w1 = 1;
        w2 = 1;
    end
    A(1,:) = (x1(1).*P1(3,1:3) - P1(1,1:3));
    A(2,:) = (x1(2).*P1(3,1:3) - P1(2,1:3));
    A(3,:) = (x2(1).*P2(3,1:3) - P2(1,1:3));
    A(4,:) = (x2(2).*P2(3,1:3) - P2(2,1:3));
    b = [P1(1,4) - P1(3,4);P1(2,4)-P1(3,4);P2(1,4)-P2(3,4);P2(2,4)-P2(3,4)];
    
    A(1:2,:) = A(1:2,:).*w1;
    A(3:4,:) = A(3:4,:).*w2;
    
    b(1:2,:) = b(1:2,:).*w1;
    b(3:4,:) = b(3:4,:).*w2;
    
    P = A\b;
end
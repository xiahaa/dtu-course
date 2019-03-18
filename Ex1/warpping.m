function Iw = warpping(I,H)
% function Iw = warpping(I,H)
%
% apply a back warpping of image I with the warpping transformation
% matrix being H. This is to say:
% U' = H*U
% where U, U' = [u,v,w]'; U' will be normalized to be a image coordinate.
% Then nearest neighboring interpolation will be used.
%
% This function has been vectorized for boosting the performance (10x faster). 
% But generally what it does is the same as:
%     %% forward 
%     for i = 1:Mn
%         for j = 1:Nn
%             new = inv(H)*[j;i;1];
%             id = floor(new./new(3))+1;
%             Iw(i,j,:) = I(id(2),id(1),:);
%         end
%     end
%
% Author: xiahaa@space.dtu.dk
% Disclaimer: This code comes with no guarantee at all and its author
%   is not liable for any damage that its utilization may cause.

    [M,N,D] = size(I);
    q = [[1;1;1] [N;1;1] [1;M;1] [N;M;1]];
    qm = H*q;
    qm = qm ./ qm(3,:);
    miny = min(qm(2,:));
    maxy = max(qm(2,:));
    minx = min(qm(1,:));
    maxx = max(qm(1,:));
    
    Mn = round(maxy - miny) + 1;
    Nn = round(maxx - minx) + 1;
    Iw = zeros(Mn, Nn, D);

    %% forward fast
%     [X,Y] = ind2sub([M,N],1:M*N);
%     NEWXY = H*[Y;X;ones(1,numel(X))];
%     NEWXY = floor(NEWXY./NEWXY(3,:))+1;
%     IND = sub2ind([Mn,Nn],NEWXY(2,:),NEWXY(1,:));
%     for i = 1:D
%         Iw(IND+(i-1)*Mn*Nn) = I([1:M*N]+(i-1)*M*N);
%     end
    
    %% backward fast
    [X,Y] = ind2sub([Mn,Nn],1:Mn*Nn);
    NEWXY = inv(H)*[Y;X;ones(1,numel(X))];
    NEWXY = floor(NEWXY./NEWXY(3,:));
    VALID = NEWXY(1,:) > 0 & NEWXY(1,:) <= N & NEWXY(2,:) > 0 & NEWXY(2,:) <= M;
    IND = 1:Mn*Nn;
    IND = IND(VALID);
    NEWIND = sub2ind([M,N],NEWXY(2,VALID),NEWXY(1,VALID));
    for i = 1:D
        Iw(IND+(i-1)*Mn*Nn) = I(NEWIND+(i-1)*M*N);
    end
end
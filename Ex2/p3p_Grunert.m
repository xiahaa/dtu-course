function [R, t] = p3p_Grunert(P, q, K)
%function [R, t] = p3p_fisch(P, q, K)
%
% Implementation of the P3P algorithm for absolute orientation proposed in
% Haralick B M, Lee C N, Ottenberg K, et al. 
% Review and analysis of solutions of the three point perspective pose estimation problem[J]. 
% International journal of computer vision, 1994, 13(3): 331-356.
%
%   Inputs:
%       P: 3D points in World coordinate system - 3xN.
%       q: 2D points in image, - 3xN (homogeneous) or 2xN (inhomogeneous).
%       K: camera intrinsics.
%   Outputs:
%       R: rotation matrix, 3x3xN, N is the number of obtained solutions.
%       t: translation vector, 3x1xN, N is the number of obtained solutions.
%
% Author: xiahaa@space.dtu.dk
% Disclaimer: This code comes with no guarantee at all and its author
%   is not liable for any damage that its utilization may cause.

    n = size(q,2);
    if size(q,1) == 2
        q = [q;ones(1,n)];% to homogeneous
    end

    %% page 332-333, find a b c
    a = norm(P(:,2)-P(:,3));
    b = norm(P(:,1)-P(:,3));
    c = norm(P(:,1)-P(:,2));
    
    %% page 333, find j1, j2, j3, alpha, beta, gamma
    % to normalized coordinates, then f equals to 1
    % other option would be
    % f = K(1,1); px = K(1,3); py = K(2,3);
    % q(3,:) = f;
    % q(1,:) = q(1,:) - px;
    % q(2,:) = q(2,:) - py;
    q = K\q;
    j1 = q(:,1) ./ norm(q(:,1));
    j2 = q(:,2) ./ norm(q(:,2));
    j3 = q(:,3) ./ norm(q(:,3));

    cos_alpha = dot(j2,j3);
    cos_beta = dot(j1,j3);
    cos_gamma = dot(j1,j2);
    
    %% some constants
    a2 = a*a;
    b2 = b*b;
    c2 = c*c;
    %% page 334, find coefficients for the 4th order polynomial
    A4 = ((a2-c2)/b2-1)^2 - 4*c2/b2*cos_alpha^2;
    A3 = 4*( (a2-c2)/b2*(1-(a2-c2)/b2)*cos_beta - (1-(a2+c2)/b2)*cos_alpha*cos_gamma + 2*c2/b2*cos_alpha^2*cos_beta );
    A2 = 2*(((a2-c2)/b2)^2 - 1 + 2*((a2-c2)/b2)^2*cos_beta^2+2*((b2-c2)/b2)*cos_alpha^2-4*((a2+c2)/b2)*cos_alpha*cos_beta*cos_gamma+2*((b2-a2)/b2)*cos_gamma^2);
    A1 = 4*(-(a2-c2)/b2*(1+(a2-c2)/b2)*cos_beta+2*a2/b2*cos_gamma^2*cos_beta-(1-(a2+c2)/b2)*cos_alpha*cos_gamma);
    A0 = (1+(a2-c2)/b2)^2-4*a2/b2*cos_gamma^2;
    
    %% eq9, the 4th order polynomial
    Poly = [A4, A3, A2, A1, A0];
    root = roots(Poly);
    
    %% find real roots
    nvalid = (abs(imag(root))<1e-6);
    v = real(root(nvalid));
    
    %% some constants
    cc1 = -1+(a2-c2)/b2;
    cc2 = -2*((a2-c2)/b2)*cos_beta;
    cc3 = 1+(a2-c2)/b2;
    
    s123 = [];
    for i = 1:length(v)
        %% eq8
        u = (cc1*v(i)^2+cc2*v(i)+cc3)/(2*(cos_gamma-v(i)*cos_alpha));
        
        %% eq5
        s12 = c2 / (1+u*u-2*u*cos_gamma);
        if s12 < 0; continue; end
        %% eq4
        s1 = sqrt(s12);
        s2 = u*s1;
        s3 = v(i)*s1;
        
        if s1 < 0 || s2 < 0 || s3 < 0; continue; end 
        
        s123 = [s123;[s1 s2 s3]];
    end
    
    %%
    R = zeros(3,3,size(s123,1));
    t = zeros(3,1,size(s123,1));
    
    for i = 1:size(s123,1)
        %% page 333, pi = si*ji 
        Q = [s123(i,1)*j1 s123(i,2)*j2 s123(i,3)*j3];
        
        %% Appendix I A Simple Linear Solution for the Absolute Orientation 
        [Ropt,topt] = lssol(P(:,1:3), Q);
        
        R(:,:,i) = Ropt;
        t(:,:,i) = topt;
    end
end

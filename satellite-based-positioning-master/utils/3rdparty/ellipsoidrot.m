function [xx,yy,zz] = ellipsoidrot(xc,yc,zc,xr,yr,zr,Q,n)
    %ELLIPSOID Generate ellipsoid.
    %
    % [X,Y,Z] = ELLIPSOID(XC,YC,ZC,XR,YR,ZR,Q,N) generates three
    % (N+1)-by-(N+1) matrices so that SURF(X,Y,Z) produces a rotated
    % ellipsoid with center (XC,YC,ZC) and radii XR, YR, ZR.
    %
    % [X,Y,Z] = ELLIPSOID(XC,YC,ZC,XR,YR,ZR,Q) uses N = 20.
    %
    % ELLIPSOID(...) and ELLIPSOID(...,N) with no output arguments
    % graph the ellipsoid as a SURFACE and do not return anything.
    %
    % The ellipsoidal data is generated using the equation after rotation with
    % orthogonal matrix Q:
    %
    % (X-XC)?2 (Y-YC)?2 (Z-ZC)?2
    % -------- + -------- + -------- = 1
    % XR?2 YR?2 ZR?2
    %
    % See also SPHERE, CYLINDER.
    % Modified by Allan Aasbjerg Nielsen (2004) after
    % Laurens Schalekamp and Damian T. Packer
    % Copyright 1984-2002 The MathWorks, Inc.
    % $Revision: 1.7 $ $Date: 2002/06/14 20:33:49 $
    error(nargchk(7,8,nargin));
    if nargin == 7
        n = 20;
    end
    [x,y,z] = sphere(n);
    x = xr*x;
    y = yr*y;
    z = zr*z;
    xvec = Q*[reshape(x,1,(n+1)^2); reshape(y,1,(n+1)^2); reshape(z,1,(n+1)^2)];
    x = reshape(xvec(1,:),n+1,n+1)+xc;
    y = reshape(xvec(2,:),n+1,n+1)+yc;
    z = reshape(xvec(3,:),n+1,n+1)+zc;
    if(nargout == 0)
        surf(x,y,z)
        % surfl(x,y,z)
        % surfc(x,y,z)
        axis equal
        %shading interp
        %colormap gray
    else
       xx = x;
       yy = y;
       zz = z;
    end
end
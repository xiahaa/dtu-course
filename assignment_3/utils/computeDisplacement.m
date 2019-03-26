function displacement = computeDisplacement(im, curve, cin, cout)
%Compute displacements for curves.
	% interpolate curve intensity
    Icurve = biInterIntensity(im, curve);
    % scalar
    s = (cin - cout).*(2*Icurve - cin - cout);
    % project to normal direction
    normals = computeNormal(curve);
    % debug; draw quiver to see if normals are computed correctly
%     normals = sign(s).*normals;
%     imshow(im);hold on;plot(curve(2,:),curve(1,:),'r-');
%     quiver(curve(2,:),curve(1,:),normals(2,:),normals(1,:));

    displacement = s.*normals;
end
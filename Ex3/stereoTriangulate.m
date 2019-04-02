
function xyz = stereoTriangulate(disp, f, cx, cy, b,id)
    z = f * b ./ disp(id);
    [yy,xx] = find(id);
    x = (xx - cx).*z / f;
    y = (yy - cy).*z / f;
    xyz = [x y z];
end

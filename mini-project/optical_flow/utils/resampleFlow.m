function [flow_u, flow_v] = resampleFlow(u, v, sz)
% instead of using hard scale, here I do the flow updating using bilinear
% interpolation and actual ratio between cascaded layers.
    ratio = sz(1) / size(u,1);
    flow_u = imresize(u,sz,'bilinear').*ratio;
    flow_v = imresize(v,sz,'bilinear').*ratio;
end
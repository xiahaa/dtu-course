function noise=perlin(values,x,y,z)
    if(numel(values)~=512)
        values=randperm(256)-1;
        values=[values values];
    end
    x=abs(x);
    X=bitand(floor(x),255);
    x=x-floor(x);
    u=fade(x);
    A=values(1+X);
    noise=linterp(u,grad1d(values(1+X),x),grad1d(values(1+X+1),x-1));
    if(nargin>2)
        y=abs(y);
        Y=bitand(floor(y),255);
        y=y-floor(y);
        v=fade(y);
        A=A+Y;
        B=values(1+X+1)+Y;
        noise=linterp(u,linterp(u,grad2d(values(1+A),x,y),grad2d(values(1+B),x-1,y)),linterp(u,grad2d(values(1+A+1),x,y-1),grad2d(values(1+B+1),x-1,y-1)));
    end
    if(nargin>3)
        z=abs(z);
        Z=bitand(floor(z),255);
        z=z-floor(z);
        w=fade(z);
        AA=values(1+A)+Z;
        AB=values(1+A+1)+Z;
        BA=values(1+B)+Z;
        BB=values(1+B+1)+Z;
        noise=linterp(  w, ... 
                        linterp(v, ... 
                                linterp(u, ... 
                                        grad3d(values(1+AA),x,y,z), ... 
                                        grad3d(values(1+BA),x-1,y,z)), ...
                                linterp(u, ...
                                        grad3d(values(1+AB),x,y-1,z), ...
                                        grad3d(values(1+BB),x-1,y-1,z))), ...
                        linterp(v, ...
                                linterp(u, ... 
                                        grad3d(values(1+AA+1),x,y,z-1), ... 
                                        grad3d(values(1+BA+1),x-1,y,z-1)), ...
                                linterp(u, ...
                                        grad3d(values(1+AB+1),x,y-1,z-1), ...
                                        grad3d(values(1+BB+1),x-1,y-1,z-1))));
    end
end
% 
%  function res = calcNoise(x,y,z,p)
%     X = bitand(floor(x),255,'uint8');
%     Y = bitand(floor(y),255,'uint8');
%     Z = bitand(floor(z),255,'uint8');
%     
%     x = x-floor(x);
%     y = y-floor(y);
%     z = z-floor(z);
% 
%     % Compute fade curves for each of x, y, z
%     u = fade(x);
%     v = fade(y);
%     w = fade(z);
%     
%     % Hash coordinates of the 8 cube corners
%     A  = p(X+1) + Y;
%     AA = p(A+1) + Z;
%     AB = p(A + 2) + Z;
%     B  = p(X + 2) + Y;
%     BA = p(B) + Z;
%     BB = p(B + 1) + Z;
%     
%     % Add blended results from 8 corners of cube
%     res = lerp(w, lerp(v, lerp(u, grad(p(AA+1), x, y, z), grad(p(BA+1), x - 1, y, z)), ...
%               lerp(u, grad(p(AB+1), x, y - 1, z), grad(p(BB+1), x - 1, y - 1, z))), ...
%                 lerp(v, lerp(u, grad(p(AA+2), x, y, z - 1), grad(p(BA + 2), x - 1, y, z - 1)), ...
%                     lerp(u, grad(p(AB + 2), x, y - 1, z - 1), grad(p(BB + 2), x - 1, y - 1, z - 1))));
%     res = (res + 1.0) / 2.0;
%  end

function l=linterp(t,a,b)
    l=a+t*(b-a);
end

function t=fade(t)
    t=6*t^5-15*t^4+10*t^3;
end

function g=grad1d(hash,x)
    if(bitand(hash,1))
        g=-x;
    else
        g=x;
    end
end

function g=grad2d(hash,x,y)
    h=bitand(hash,3);
    if(bitand(h,2))
        u=-x;
    else
        u=x;
    end
    if(bitand(h,1))
        v=-y;
    else
        v=y;
    end
    g=u+v;
end

function g=grad3d(hash,x,y,z)
    h=bitand(hash,15);
    if(h<8)
        u=x;
    else
        u=y;
    end
    if(h<4)
        v=y;
    elseif(h==12 || h==14)
        v=x;
    else
        v=z;
    end
    if(bitand(h,1))
        if(bitand(h,2))
            g=-u-v;
        else
            g=-u+v;
        end
    else
        if(bitand(h,2))
            g=u-v;
        else
            g=u+v;
        end
    end
end
% 
% function n = grad(hash,x,y,z)
%     h = bitand(hash,15,'int32');%    hash & 15;
%     if h < 8
%         u = x;
%     else
%         u = y;
%     end
%     if h < 4
%         v = y;
%     else
%         if h == 12 || h == 14
%             v = x;
%         else
%             v = z;
%         end
%     end
%     if bitand(h,1,'int32') == 0
%         if bitand(h,2,'int32') == 0
%             n = u + v;
%         else
%             n = u - v;
%         end
%     else
%         if bitand(h,2,'int32') == 0
%             n = -u + v;
%         else
%             n = -u - v;
%         end
%     end
% end
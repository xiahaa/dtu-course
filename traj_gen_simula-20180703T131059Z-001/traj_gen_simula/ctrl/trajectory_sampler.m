function [ desired_state ] = trajectory_sampler(t, qn, current_state, trajCoeff, trajTime, mode)
% TRAJECTORY_SAMPLER: sample a desired state from the generated trajectory
%
% desired_state: Contains all the information that is passed to the
% controller, as in phase 2
%
% It is suggested to use "persistent" variables to store map and path
% during the initialization call of trajectory_generator, e.g.
% persistent map0 path0
% map0 = map;
% path0 = path;
persistent total_time X ts sampleMode;
if numel(t) == 0 | numel(qn) == 0
   total_time = trajTime(end);
   ts = trajTime;
   X = trajCoeff;
   sampleMode = mode;
   return
end

numCoeff = 6;

if sampleMode == 1 || nargin == 2
    if t >= total_time
        k = numel(ts)-1;
        scalar = 1./(ts(end)-ts(end-1));
        tnorm = 1;
    else
        % 5th order minimum snap trajectory
    %     t
        k = find(ts<=t);
        k = k(end);
        scalar = 1./(ts(k+1)-ts(k));
        tnorm = (t - ts(k))*scalar;
    %     scalar
    %     tnorm
    end
elseif sampleMode == 2
    P = [91.4417   -0.0000   -0.0000    9.5495   -0.0000   -0.0000
        -0.0000  128.6153   -0.0000   -0.0000   10.3090   -0.0000
        -0.0000   -0.0000  111.9488   -0.0000   -0.0000    6.8348
         9.5495   -0.0000   -0.0000   12.2950   -0.0000   -0.0000
        -0.0000   10.3090   -0.0000   -0.0000   13.1498   -0.0000
        -0.0000   -0.0000    6.8348   -0.0000   -0.0000    6.3511];
    
    persistent coarseSamplesP coarseSamplesT
    if numel(coarseSamplesP) == 0 || numel(coarseSamplesT) == 0
        tss = ts(1):0.01:ts(end);tss = tss';
        for j = 1:numel(tss)
            k = find(ts<=tss(j),1,'last');
            k(k==numel(ts)) = k(k==numel(ts)) - 1;
            scalar = 1./(ts(k+1));
            tnorm = (tss(j) - ts(k))*scalar;
            
            polyx = X((k-1)*numCoeff+1:k*numCoeff,1);
            polyy = X((k-1)*numCoeff+1:k*numCoeff,2);
            polyz = X((k-1)*numCoeff+1:k*numCoeff,3);
            
            xsamples = polyx(1) + polyx(2).*tnorm + polyx(3).*tnorm.^2 + polyx(4).*tnorm.^3 + polyx(5).*tnorm.^4 + polyx(6).*tnorm.^5;
            ysamples = polyy(1) + polyy(2).*tnorm + polyy(3).*tnorm.^2 + polyy(4).*tnorm.^3 + polyy(5).*tnorm.^4 + polyy(6).*tnorm.^5;
            zsamples = polyz(1) + polyz(2).*tnorm + polyz(3).*tnorm.^2 + polyz(4).*tnorm.^3 + polyz(5).*tnorm.^4 + polyz(6).*tnorm.^5;
            
            vxsamples = polyx(2) + polyx(3).*tnorm.*2 + polyx(4).*3.*tnorm.^2 + polyx(5).*4.*tnorm.^3 + polyx(6).*5.*tnorm.^4;
            vysamples = polyy(2) + polyy(3).*tnorm.*2 + polyy(4).*3.*tnorm.^2 + polyy(5).*4.*tnorm.^3 + polyy(6).*5.*tnorm.^4;
            vzsamples = polyz(2) + polyz(3).*tnorm.*2 + polyz(4).*3.*tnorm.^2 + polyz(5).*4.*tnorm.^3 + polyz(6).*5.*tnorm.^4;
            
            coarseSamplesP(j,:) = [xsamples ysamples zsamples vxsamples vysamples vzsamples];
            coarseSamplesT(j) = tss(j);
        end
    end
    
    cs = current_state(1:6)';
    
    err = (repmat(cs, size(coarseSamplesP,1),1) - coarseSamplesP);
    dist = diag(err * P * err');
    [~,minid] = min(dist);
    
    tv = coarseSamplesT(minid);
    if tv >= total_time
        k = numel(ts)-1;
        scalar = 1./(ts(end)-ts(end-1));
        tnorm = 1;
    else
        % 5th order minimum snap trajectory
        k = find(ts<=t);
        k = k(end);
        scalar = 1./(ts(k+1)-ts(k));
        tnorm = (tv - ts(k))*scalar;
    end 
elseif sampleMode == 3
    
end

polyx = X((k-1)*numCoeff+1:k*numCoeff,1);
polyy = X((k-1)*numCoeff+1:k*numCoeff,2);
polyz = X((k-1)*numCoeff+1:k*numCoeff,3);

xsamples = polyx(1) + polyx(2).*tnorm + polyx(3).*tnorm.^2 + polyx(4).*tnorm.^3 + polyx(5).*tnorm.^4 + polyx(6).*tnorm.^5;
ysamples = polyy(1) + polyy(2).*tnorm + polyy(3).*tnorm.^2 + polyy(4).*tnorm.^3 + polyy(5).*tnorm.^4 + polyy(6).*tnorm.^5;
zsamples = polyz(1) + polyz(2).*tnorm + polyz(3).*tnorm.^2 + polyz(4).*tnorm.^3 + polyz(5).*tnorm.^4 + polyz(6).*tnorm.^5;

vxsamples = polyx(2) + polyx(3).*tnorm.*2 + polyx(4).*3.*tnorm.^2 + polyx(5).*4.*tnorm.^3 + polyx(6).*5.*tnorm.^4;
vysamples = polyy(2) + polyy(3).*tnorm.*2 + polyy(4).*3.*tnorm.^2 + polyy(5).*4.*tnorm.^3 + polyy(6).*5.*tnorm.^4;
vzsamples = polyz(2) + polyz(3).*tnorm.*2 + polyz(4).*3.*tnorm.^2 + polyz(5).*4.*tnorm.^3 + polyz(6).*5.*tnorm.^4;

axsamples = polyx(3).*2 + polyx(4).*6.*tnorm.^1 + polyx(5).*12.*tnorm.^2 + polyx(6).*20.*tnorm.^3;
aysamples = polyy(3).*2 + polyy(4).*6.*tnorm.^1 + polyy(5).*12.*tnorm.^2 + polyy(6).*20.*tnorm.^3;
azsamples = polyz(3).*2 + polyz(4).*6.*tnorm.^1 + polyz(5).*12.*tnorm.^2 + polyz(6).*20.*tnorm.^3;

vxsamples = vxsamples.*scalar;
vysamples = vysamples.*scalar;
vzsamples = vzsamples.*scalar;

axsamples = axsamples.*(scalar^2);
aysamples = aysamples.*(scalar^2);
azsamples = azsamples.*(scalar^2);

pos = [xsamples ysamples zsamples];
vel = [vxsamples vysamples vzsamples];
acc = [axsamples aysamples azsamples];

yaw = 0;
yawdot = 0;

% =================== Your code ends here ===================

desired_state.pos = pos(:);
desired_state.vel = vel(:);
desired_state.acc = acc(:);
desired_state.yaw = yaw;
desired_state.yawdot = yawdot;

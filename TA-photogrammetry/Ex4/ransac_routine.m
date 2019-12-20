
function result = ransac_routine(x, ransac)
% function result = ransac_routine(x, ransac)
%
% General RANSAC routine.
%   Inputs:
%       x: data.
%  ransac: RANSAC options.
%        pinlier: estimated inlier percentage.
%        estt_fun: function handler to model fitting.
%        eval_fun: function handler to model evaluation.
%        maxiter: maximum iteration number.
%        threshold: inlier threshold.
%        inliers: estimated inliers.
%        minimumset: minimum number of samples for model fitting.
%   Outputs:
%   result: RANSAC results.
%        params: estimated parameters.
%       inliers: estimated inliers.
%
% Author: xiahaa@space.dtu.dk
% Disclaimer: This code comes with no guarantee at all and its author
%   is not liable for any damage that its utilization may cause.
    iter = 1;
    n = size(x,2);
    while iter < ransac.maxiter
        % sample
        id = randperm(n, ransac.minimumset);
        % estimate
        params = ransac.estt_fun(x(:,id));
        % consensus
        errs = ransac.eval_fun(x, params);
        % verify
        inliers = errs < ransac.threshold;
        num_inlier = sum(inliers);
        % update
        if num_inlier > sum(ransac.inliers)
            ransac.inliers = inliers;
            pin = num_inlier / n;
            ransac.maxiter = round(log(1-ransac.pinlier)/log(1-pin^(ransac.minimumset)));
        end
        iter = iter + 1;
    end
    % refinement
    result.params = ransac.estt_fun(x(:,ransac.inliers));
    result.inliers = ransac.inliers;
end
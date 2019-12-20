function F = calcF(t, s, qn, controlhandle, trajhandle, params)
%
% INPUTS:
% t             - 1 x 1, time
% s             - 13 x 1, state vector = [x, y, z, xd, yd, zd, qw, qx, qy, qz, p, q, r]
% qn            - quad number (used for multi-robot simulations)
% controlhandle - function handle of your controller
% trajhandle    - function handle of your trajectory generator
% params        - struct, output from crazyflie() and whatever parameters you want to pass in
%
% OUTPUTS:
% F             - 1 x 1, thrust force
%
% NOTE: You should not modify this function
% See Also: quadEOM_readonly, crazyflie

% convert state to quad stuct for control
qd{qn} = stateToQd(s);

% Get desired_state
desired_state = trajhandle(t, qn);

% The desired_state is set in the trajectory generator
qd{qn}.pos_des      = desired_state.pos;
qd{qn}.vel_des      = desired_state.vel;
qd{qn}.acc_des      = desired_state.acc;
qd{qn}.yaw_des      = desired_state.yaw;
qd{qn}.yawdot_des   = desired_state.yawdot;

% get control outputs
[F, M, trpy, drpy] = controlhandle(qd, t, qn, params);

% compute derivative
sdot = quadEOM_readonly(t, s, F, M, params);
F = sdot(end);
end

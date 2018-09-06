function ctrl_entry(h_fig, trajCoeff, trajTime)
    addpath('ctrl/utils') % todo
    addpath('ctrl/trajectories')% ok

    % You can change trajectory here

    % trajectory generator
    % trajhandle = @step;
    % trajhandle = @circle;
    trajhandle = @trajectory_sampler;

    % controller
    controlhandle = @controller;

    % real-time 
    real_time = true;

    % *********** YOU SHOULDN'T NEED TO CHANGE ANYTHING BELOW **********
    % number of quadrotors
    nquad = 1;

    % max time
    time_tol = 25;

    % parameters for simulation
    params = crazyflie();% todo, assign maximum thrust and angle

    trajectory_sampler([], [], trajCoeff, trajTime);
    findVirtualTime([],[], trajCoeff, trajTime, params.sampleType);
    
%% **************************** FIGURES *****************************
% fprintf('Initializing figures...\n')
% h_fig = figure;
h_3d = gca;
axis equal
grid on
view(3);
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]')
    quadcolors = lines(nquad);
set(gcf,'Renderer','OpenGL')
% set(gcf, 'Position', [1200, 1200, 800, 800]);

    %% *********************** INITIAL CONDITIONS ***********************
    fprintf('Setting initial conditions...\n')
    max_iter  = 5000;      % max iteration
    starttime = 0;         % start of simulation in seconds
    tstep     = 0.01;      % this determines the time step at which the solution is given
    cstep     = 0.05;      % image capture time interval
    nstep     = cstep/tstep;
    time      = starttime; % current time
    err = []; % runtime errors
    for qn = 1:nquad
        % Get start and stop position
        des_start = trajhandle(0, qn);
        des_stop  = trajhandle(inf, qn);
        stop{qn}  = des_stop.pos;
        x0{qn}    = init_state( des_start.pos, 0 );
        xtraj{qn} = zeros(max_iter*nstep, length(x0{qn}));
        ttraj{qn} = zeros(max_iter*nstep, 1);% storage mem for internal states
    end

    x         = x0;        % state

    pos_tol   = 0.01;
    vel_tol   = 0.01;

    %% ************************* RUN SIMULATION *************************
    OUTPUT_TO_VIDEO = 0;
    if OUTPUT_TO_VIDEO == 1
        v = VideoWriter('perlin.avi');
        open(v)
    end

    fprintf('Simulation Running....')
    % Main loop
    for iter = 1:max_iter
        iter;
        timeint = time:tstep:time+cstep;
        tic;
        % Iterate over each quad
        for qn = 1:nquad
            % Initialize quad plot
            if iter == 1
                QP{qn} = QuadPlot(qn, x0{qn}, 0.1, 0.04, quadcolors(qn,:), max_iter, h_3d);
                desired_state = trajhandle(time, qn);
                QP{qn}.UpdateQuadPlot(x{qn}, [desired_state.pos; desired_state.vel], time);
                h_title = title(sprintf('iteration: %d, time: %4.2f', iter, time));
            end

            % Run simulation
            tic
            [tsave, xsave] = ode45(@(t,s) quadEOM(t, s, qn, controlhandle, trajhandle, params), timeint, x{qn});
            toc
            disp('ode ending!');
            F = calcF(timeint(end), x{qn}, qn, controlhandle, trajhandle, params);
            xsave(end,14) = F;

            if time > params.holdTime
                x{qn}    = xsave(end, :)';
            end

            % Save to traj
            xtraj{qn}((iter-1)*nstep+1:iter*nstep,:) = xsave(1:end-1,:);
            ttraj{qn}((iter-1)*nstep+1:iter*nstep) = tsave(1:end-1);

            % Update quad plot
            desired_state = trajhandle(time + cstep, qn);
            QP{qn}.UpdateQuadPlot(x{qn}, [desired_state.pos; desired_state.vel], time + cstep);
            set(h_title, 'String', sprintf('iteration: %d, time: %4.2f', iter, time + cstep))
            if OUTPUT_TO_VIDEO == 1
                im = frame2im(getframe(gcf));
                writeVideo(v,im);
            end
            
            if params.sampleType == 1 
                time = time + cstep; % Update simulation time
            elseif params.sampleType == 2
                time = findVirtualTime([], x{qn}(1:6)) + cstep;
            else
                desired_state1 = trajhandle(time, qn);
                time = findVirtualTime([], x{qn}(1:6));
            end
            
        end
        t = toc;
        % Check to make sure ode45 is not timing out
        if(t> cstep*500)
            err = 'Ode45 Unstable';
            break;
        end

        % Pause to make real-time
        if real_time && (t < cstep)
            pause(cstep - t);
        end

        % Check termination criteria
        if terminate_check(x, time, stop, pos_tol, vel_tol, time_tol)
            break
        end
    end

    if OUTPUT_TO_VIDEO == 1
        close(v);
    end

    %% ************************* POST PROCESSING *************************
    % Truncate xtraj and ttraj
    for qn = 1:nquad
        xtraj{qn} = xtraj{qn}(1:iter*nstep,:);
        ttraj{qn} = ttraj{qn}(1:iter*nstep);
    end

    % Plot the saved position and velocity of each robot
    for qn = 1:nquad
        % Truncate saved variables
        QP{qn}.TruncateHist();
        % Plot position for each quad
        h_pos{qn} = figure('Name', ['Quad ' num2str(qn) ' : position']);
        plot_state(h_pos{qn}, QP{qn}.state_hist(1:3,:), QP{qn}.time_hist, 'pos', 'vic');
        plot_state(h_pos{qn}, QP{qn}.state_des_hist(1:3,:), QP{qn}.time_hist, 'pos', 'des');
        % Plot velocity for each quad
        h_vel{qn} = figure('Name', ['Quad ' num2str(qn) ' : velocity']);
        plot_state(h_vel{qn}, QP{qn}.state_hist(4:6,:), QP{qn}.time_hist, 'vel', 'vic');
        plot_state(h_vel{qn}, QP{qn}.state_des_hist(4:6,:), QP{qn}.time_hist, 'vel', 'des');

        qr = QP{qn}.state_hist(7:10,:);
    %     qd = QP{qn}.state_des_hist(7:10,:);
        for j = 1:size(qr,2)
            Rr = QuatToRot(qr(:,j));
    %         Rd = QuatToRot(qd(:,j));
            [rollr,pitchr,yawr] = RotToRPY_ZXY(Rr);
    %         [rolld,pitchd,yawd] = RotToRPY_ZXY(Rd);
            rrs(j) = rollr * 180/pi;
            prs(j) = pitchr * 180/pi;
    %         rds(j) = rolld;
    %         pds(j) = pitchd;
        end

        % Plot roll for each quad
        r{qn} = figure('Name', ['Quad ' num2str(qn) ' : roll']);
        plot_state(r{qn}, rrs, QP{qn}.time_hist, 'roll', 'vic');
    %     plot_state(r{qn}, rds, QP{qn}.time_hist, 'roll', 'des');
        % Plot pitch for each quad
        p{qn} = figure('Name', ['Quad ' num2str(qn) ' : pitch']);
        plot_state(p{qn}, prs, QP{qn}.time_hist, 'pitch', 'vic');
    %     plot_state(p{qn}, pds, QP{qn}.time_hist, 'pitch', 'des');
        % Plot thrust for each quad
        thrust{qn} = figure('Name', ['Quad ' num2str(qn) ' : thrust']);
        plot_state(thrust{qn}, QP{qn}.state_hist(14,:), QP{qn}.time_hist, 'thrust', 'vic');
    %     plot_state(thrust{qn}, QP{qn}.state_des_hist(14,:), QP{qn}.time_hist, 'thrust', 'des');
    end
    if(~isempty(err))
        error(err);
    end

    fprintf('finished.\n')

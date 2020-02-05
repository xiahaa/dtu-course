%% I Forward model and model parameters
%  =====================================

% 1. Defining the dynamic model for generating the data
% =====================================================
% 
% Air density model
gamma = 1.754;
eta   = 1.39e-4;
% Vectorize: Write the functions so that they allow an input of several columns
rho = @(h) gamma*exp(-eta*h);
g = 9.81;
% x = [h,v,beta]
TrackingDynamics = @(t,x) [-x(2,:);g*(1-rho(x(1,:)).*(x(2,:).^2)./(2*x(3,:)));zeros(1,size(x,2))];

% 2. Generate the data
% ====================
% Initial values
h_0       = 61000;
v_0       = 3048;
beta_true = 19161;
x0        = [h_0;v_0;beta_true];

% Time step, time span
dt     = 0.1;
t_span = dt*(0:300);
nt     = length(t_span);

[t,x] = ode45(TrackingDynamics,t_span,x0);

% Noise model
sigma = 500;
b = x(:,1) + sigma*randn(nt,1);
% Plot the data
figure(1)
plot(t,b,'k-','LineWidth',2)
set(gca,'FontSize',15)
xlabel('time [s]','FontSize',15)
ylabel('altitude [m]','FontSize',15)


% Variances of the innovations
sigma_h    = 100;
sigma_v    = 100;
sigma_beta = 5;

D_iv = diag([sigma_h,sigma_v,sigma_beta]);  % Square root of the innovation matrix

%% II Particle filtering
%  ======================

% 1. Generate the initial particle cloud
% ======================================

h_init    = b(1);
v_init    = 3000;
beta_init = 20000;
x_init    = [h_init;v_init;beta_init];  

std_h    = 500;
std_v    = 200;
std_beta = 2500;
STD      = [std_h;std_v;std_beta];

N = 5000; % Number of particles

% Current particle cloud and the weights

X_c = x_init + diag(STD)*randn(3,N);
W_c = 1/N*ones(1,N);


% 2. Particle filtering algorithm
% ===============================

x_mean      = zeros(3,nt); % Estimated mean values 
x_std       = zeros(3,nt); % Estimated standard deviations
x_mean(:,1) = x_init;
x_std(:,1)  = STD;

for j = 1:nt-1
    % 2 a. Propagate the current cloud using forward Euler
    % ====================================================
    X_pred = X_c + dt*TrackingDynamics(0,X_c);
   
    % 2 b. Compute the fitness weights
    % ================================
    W_pred = W_c.*exp(-1/(2*sigma^2)*(X_pred(1,:) - b(j+1)).^2);
    W_pred = (1/sum(W_pred))*W_pred;
    % 2 c. Draw with replacement using the fitness weights
    % ====================================================
    [W_sort,I_sort] = sort(W_pred,'descend');
    I_draws = zeros(1,N);
    for ell = 1:N
        k   = 0;
        Phi = 0;
        xi  = rand;
        while Phi < xi
            k = k + 1;
            Phi = Phi + W_sort(k);
        end
        I_draws(ell) = I_sort(k);
    end
    X_pred = X_pred(:,I_draws);
    
    % 2 d. Add innovation
    % ===================
    X_c = X_pred + D_iv*randn(3,N);
    % 2 e. Update the weights
    % =======================
    W_c = exp(-1/(2*sigma^2)*((X_c(1,:) - b(j+1)).^2 - (X_pred(1,:) - b(j+1)).^2));
    W_c = (1/sum(W_c))*W_c;
    % 2 f. Estimate the particle mean and standard deviation
    % ======================================================
    x_mean(:,j+1) = X_c*W_c';
    X_centered = X_c - x_mean(:,j+1)*ones(1,N);
    CC = X_centered*diag(W_c)*X_centered';
    x_std(:,j+1) = sqrt(diag(CC));
end


% Plot the one standard deviation envelope of each component, indicating the data in the
% first one. The true trajectories are is in red
figure(2)
h_low  = x_mean(1,:) - x_std(1,:);
h_high = x_mean(1,:) + x_std(1,:);
fill([t_span,t_span(end:-1:1)],[h_low,h_high(end:-1:1)],[0.95,0.95,1]);
hold on
plot(t_span,x(:,1),'r-','LineWidth',2)
plot(t_span,x_mean(1,:),'k-','LineWidth',2)
plot(t_span,b,'b.','MarkerSize',10)
set(gca,'FontSize',15)
hold off

figure(3)
v_low  = x_mean(2,:) - x_std(2,:);
v_high = x_mean(2,:) + x_std(2,:);
fill([t_span,t_span(end:-1:1)],[v_low,v_high(end:-1:1)],[0.95,0.95,1]);
hold on
plot(t_span,x(:,2),'r-','LineWidth',2)
hold on
plot(t_span,x_mean(2,:),'k-','LineWidth',2)
set(gca,'FontSize',15)
hold off

figure(4)
beta_low  = x_mean(3,:) - x_std(3,:);
beta_high = x_mean(3,:) + x_std(3,:);
fill([t_span,t_span(end:-1:1)],[beta_low,beta_high(end:-1:1)],[0.95,0.95,1]);
hold on
plot(t_span,x(:,3),'r-','LineWidth',2)
hold on
plot(t_span,x_mean(3,:),'k-','LineWidth',2)
set(gca,'FontSize',15)
hold off


    
    
    
    




clc; clear; close all;
fprintf('Blending Process Simulation\n')
% Install Control System Toolbox, Symbolic Math Toolbox


%% --- Display Manipulated, Controlled, and Disturbance Variables ---
fprintf('\nManipulated Variable  : w2 \n')
fprintf('Controlled Variable   : x  \n')
fprintf('Disturbance Variable  : x1 \n')

% V is the internal condition of system which is given in starting line.
fprintf('State Variable        : V  \n\n')


%% We define function

global w1 w2 x1 x2 rho

w1 = 500;
w2 = 200;
x1 = 0.4;
x2 = 1;
rho = 900;
V = 1;

sampling_time = 1;
N = 200;

function dV = dV_dt(V,w2_current)    % In this V, w2 are variables
    global w1 rho
    dV = (w1 + w2_current - ( w1 + w2_current)) / rho;
end

function dVx = dVx_dt(V,x,w2_current)     % In this V, w2, x are variables
    global w1 x1 x2 rho
    dVx = ( w1*x1 + w2_current*x2 - ( w1 + w2_current)*x ) / ( rho*V );

end


%% --- Steady-State Simulation ---

sampling_time = 1;
N = 200;

% Time is changing from 0 to 200 with differance of sampling_time (1)
% we are starting from 0 because there is no time delay
tspan = 0:sampling_time:N;

V0 = 1;                % initial volume
x0 = 0;                % I am assuming initial fraction is zero
y0 = [V0; x0];


% Define ODE function with variables of t, y (where y have two variables
% y(1) is V and y(2 is x)) at the value w2 = 200
ode_steady = @(t, y) [dV_dt(y(1), w2); dVx_dt(y(1), y(2), w2)];


% Run ODE45 simulation over time 
[t_steady, y_steady] = ode45(ode_steady, tspan, y0);


% First cloumn y is volume and second is molel fraction
% here y(1) is v cloume and it is constant with the value of 1
V = y_steady(:,1);
x = y_steady(:,2);

figure;
plot(t_steady, x, 'y', 'LineWidth', 2);
xlabel('Time (min)');
ylabel('Mass fraction of A (x)');
title('Steady-State Identification');
grid on;



%% Open Loop Dynamic Simulation


w2_ss = w2; % store steady-state value

% we are changing value of w2 with some step change
step_changes = [0.10, 0.20];  % +10% and +20% step increases

t_step = 100;                 % time of step

out_profiles = {};

% we run a for loop for two step change
for i = 1:length(step_changes)
    
    w2_step = w2_ss * (1 + step_changes(i)); % new w2 value
    
    % Define w2(t) its means that before t=100 value of w2 remians same but
    % after t=100 we are doing some step change in w2
    w2_func = @(t) (t < t_step)*w2_ss + (t >= t_step)*w2_step;
    

    % Define ODE at w2_func which is change according to time
    ode_open_loop = @(t, y) [dV_dt(y(1), w2_func(t)); dVx_dt(y(1), y(2), w2_func(t)) ];
    
    % We simulate ODE45
    [t, y] = ode45(ode_open_loop, tspan, y0);
    
    V = y(:,1);
    x = y(:,2);
    
    % Add 1% measurement noise
    noise_var = 0.01 * mean(abs(x));
    x_noisy = x + sqrt(noise_var) * randn(size(x));
    
    % Now we are storing our all values like time, step change fraction, step change noise 
    out_profiles{i} = struct('t', t, 'x', x, 'x_noisy', x_noisy);
end

%Plot results
figure; hold on;
noisy = ['b','r'];
ideal = ['y','g'];

for i = 1:length(step_changes)
    plot(out_profiles{i}.t, out_profiles{i}.x_noisy, noisy(i), 'LineWidth', 2); hold on;
    plot(out_profiles{i}.t, out_profiles{i}.x, ideal(i), 'LineWidth', 2);
end

xlabel('Time (min)');
ylabel('Mass fraction of A in outlet (x)');
title('Open Loop Dynamic Simulation');
legend('+10% Step (noisy)', 'Ideal +10%', '+20% Step (noisy)', 'Ideal +20%');
grid on;


%% Calculate Kp and Tau by theoritacal without Step chnage


% Analytical steady-state at operating point
x_ss_theori = (w1*x1 + w2_ss*x2) / (w1 + w2_ss);

% Linearized (infinitesimal) gain and analytic tau
Kp_theori = (x2 - x_ss_theori) / (w1 + w2_ss);   % = w1*(x2-x1)/(w1+w2)^2
tau_theori = (rho * V0) / (w1 + w2_ss);

fprintf('ANALYTIC results without Step change:\n');
fprintf('  x_ss (Steady) = %.4f\n', x_ss_theori);
fprintf('  Kp            = %.6f\n', Kp_theori);
fprintf('  tau           = %.6f min\n\n', tau_theori);




%% --- Calculate Kp and Tau ---


for idx = 1:length(step_changes)
    t = out_profiles{idx}.t;
    x = out_profiles{idx}.x;

    % Steady values
    x0 = mean(x(30:99));      % We took steady state value as
    x_ss = mean(x(130:end));  % final steady
    delta_x = x_ss - x0;
    
    delta_w2 = w2_ss * step_changes(idx);
    
    % Process gain
    Kp(idx) = delta_x / delta_w2;
    
    % Time constant (63.2% rise)
    % after step chage how much time it take to become steady
    % So we find steady value after step change and substract it from at step change
    x_target = x0 + 0.632 * delta_x;
    [~ , idx_tau] = min(abs(x - x_target));
    tau(idx) = t(idx_tau) - t(t_step);


    fprintf('Step +%.0f%%:\n  x0_mean = %.6f\n  xf_mean = %.6f\n  delta_x = %.6f\n  delta_w2 = %.3f kg/min\n  Kp_step = %.6f\n  tau_step = %.4f min\n\n', ...
        step_changes(idx)*100, x0, x_ss, delta_x, delta_w2, Kp(idx), tau(idx));

end


%% Linearization and Simulation of Linear State Space Model and Stability Analysis

% Steady-state calculations
w = w1 + w2;               % steady-state outlet flow
x0 = (w1*x1 + w2*x2)/w;    % steady-state outlet composition


% Nonlinear model:
% f = dx/dt = (1/(rho*V))*(w1*(x1 - x) + w2*(x2 - x))
% here x, w2, x1 are variables. We apply Taylor's series expansion
% Linearization around steady state
% A = ∂(f)/∂x
% B1 = ∂(f)/∂w2  (manipulated variable)
% B2 = ∂(f)/∂x1  (disturbance variable)

A  = -(w1 + w2)/(rho * V0);
V1 = (x2 - x0)/(rho * V0);
V2 = w1/(rho * V0);
C  = 1;
D  = [0 0];

fprintf('---- Linearized Model Matrices ----\n');
fprintf('\n[A: %.4f , B1: %.4f, B2: %.4f]\n\n',A, V1, V2)

syms x w2 x1 real
% Nonlinear function definition
f = (1/(rho*V0)) * (w1*(x1 - x) + w2*(x2 - x));   % dx/dt expression

% Jacobians
A1  = simplify(diff(f, x));   % ∂f/∂x
V11 = simplify(diff(f, w2));  % ∂f/∂w2  -> Manipulated variable
V21 = simplify(diff(f, x1));  % ∂f/∂x1  -> Disturbance variable

% Substitute steady-state values
A_val  = double(subs(A1, [x, w2, x1], [x0, w2_ss, x1]));
V1_val = double(subs(V11, [x, w2, x1], [x0, w2_ss, x1]));
V2_val = double(subs(V21, [x, w2, x1], [x0, w2_ss, x1]));

% Display results
fprintf('---- Linearized Model (Jacobian Method) ----\n');
fprintf('A  = ∂f/∂x  = %.6f\n', A_val);
fprintf('V1 = ∂f/∂w2 = %.6f\n', V1_val);
fprintf('V2 = ∂f/∂x1 = %.6f\n', V2_val);

%%
fprintf('\n---- MISO ----\n')
% build MISO
B_val = [V1_val V2_val];  % [input: w2 , disturbance: x1]
sys_lin = ss(A_val, B_val, C, D);

% Convert to Transfer Function (each input separately)
[num1, den] = ss2tf(A_val, B_val(:,1), C, D(:,1));   % From w2 → x
[num2, ~]   = ss2tf(A_val, B_val(:,2), C, D(:,2));   % From x1 → x

G_w2 = tf(num1, den);   % Manipulated variable (input)
G_x1 = tf(num2, den);   % Disturbance variable

fprintf('\n---- Transfer Functions ----\n');
disp('From w2 → x :')
G_w2
disp('From x1 → x :')
G_x1


% Stability analysis
poles = pole(G_w2);
fprintf('\n---- Poles ----')
fprintf('\nSystem Poles: %.4f\n',poles);

if all(real(poles) < 0)
    disp('System is Stable (all poles negative)');
else
    disp('System is Unstable (pole(s) positive or zero)');
end




%% Step Response Comparison for w2 step changes (10% & 20%)
t = 0:1:200;           % time (s)
dt = 1;                % integration step
V = V0;

clear x w2 x1

w2 = 200;
x1 = 0.4;

% Define step magnitudes (in %)
colors = ['b', 'y'];   % blue for 10%, red for 20%

figure; hold on;

for i = 1:length(step_changes)

    % Nonlinear model
    x_nonlinear = x0;
    x_store = zeros(size(t));

    % we start a for loop w.r.t time
    for k = 1:length(t)
        if t(k) > 100   % there is step change
            w2_input = (1 + step_changes(i)) * w2;
        else
            w2_input = w2;
        end
        dx_dt = (1/(rho*V)) * (w1*(x1 - x_nonlinear) + w2_input*(x2 - x_nonlinear));
        x_nonlinear = x_nonlinear + dx_dt*dt;
        x_store(k) = x_nonlinear;
    end

    % Linear model 
    u = zeros(length(t),2);     % we are take two column for w2 and x1
    u(t >= 100,1) = step_changes(i) * w2;   % step change in w2
    u(:,2) = 0;                     % no change in x1
    % we apply linear simulation function. we linerize wrt to w2 and t
    [x_lin, t_lin] = lsim(sys_lin, u, t);

    % plot
    plot(t, x_store, colors(i), 'LineWidth', 1.5);
    plot(t_lin, x0 + x_lin(:,1), [colors(i) '--'], 'LineWidth', 1.5);
end

xlabel('Time (s)');
ylabel('Mass fraction of A, x');
title('Comparison: Nonlinear vs Linear Response for 10% and 20% Step in w_2');
legend('Nonlinear (10%)','Linear (10%)','Nonlinear (20%)','Linear (20%)');
grid on;

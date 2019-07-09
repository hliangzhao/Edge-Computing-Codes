%% This script simulates the Lyapunov Optimization-based Dynamic Computation Offloading (LODCO) algorithm.
%  author: Hailiang Zhao
clc, clear
opt = optimset('Display', 'none');

%% basic parameter settings (had better not change those paras)
k = 1e-28;                % effective switched capacitance (a constant decided by the chip architecture)
tau = 0.002;              % the length of time slot (in second)
phi = 0.002;              % the cost of task dropping (in second)
omega = 1e6;              % the bandwidth of MEC server (in Hz)
sigma = 1e-13;            % the noise power of the receiver (in W)
p_tx_max = 1;             % the maximum transmit power of mobile device (in W)
f_max = 1.5e9;            % the maximum CPU-cycle frequency of mobile device (in Hz)
E_max = 0.002;            % the maximum amout of battery output energy (in J)
L = 1000;                 % the input size of the computation task (in bit)
X = 737.5;                % the number of CPU cycles needed on processing one bit of task
W = L * X;                % the number of CPU cycles needed on processing one task
E_H_max = 48e-6;          % the upper bound of the energy arrive at the mobile device (in J)
p_H = E_H_max / (2*tau);  % the average Energy Harvesting (EH) power (in W)
g0 = power(10, -4);       % the path-loss constant

%% parameter control
T = 50000;                % the number of time slot (a.k.a. the size of the time horizon)
tau_d = 0.002;            % execution deadline (in second)
d = 50;                   % the distance between the mobile device and the MEC server (in meter)
E_min = 0.02e-3;          % the minimum amout of battery output energy (in J)
V = 1e-5;                 % the weight of penalty (the control parameter introduced by Lyapunov Optimization)
rho = 0.6;                % the probability that the computation task is requested

% the lower bound of perturbation parameter
E_max_hat = min(max(k * W * (f_max)^2, p_tx_max * tau), E_max);
theta = E_max_hat + V * phi / E_min;

%% allocate storage for valuable results
B = zeros(T, 1);          % the battery energy level (in J)
B_hat = zeros(T, 1);      % the virtual battery energy level ($B_hat = B - theta$)
e = zeros(T, 1);          % the amout of the harvested and stored energy (in J)
chosen_mode = zeros(T, 1);% {1: local, 2: remote, 3: drop, 4: no task request}
f = zeros(T, 1);          % the CPU-cycle frequency of local execution (in Hz)
p = zeros(T, 1);          % the transmit power of computation offloading (in W)
cost = zeros(T, 3);       % execution delay for mobile execution, MEC server execution and final choice, respectively (in second)
E = zeros(T, 3);          % energy consumption for mobile execution, MEC server execution and final choice, respectively (in J)

%% simulation begin
t = 1;
while t <= T
    disp(['===> Time slot #', num2str(t), ' <==='])
    
    %% initialization
    % generate the task request
    zeta = binornd(1, rho);
    % generate the virtual battery energy level
    B_hat(t) = B(t) - theta;
    
    %% step 1: get the optimal energy harvesting no matter whether task is requested
    E_H_t = unifrnd(0, E_H_max);
    if B_hat(t) <= 0
        e(t) = E_H_t;
    end
    
    %% step 2: get the optimal computation offloading strategy (I_m, I_s, I_d, f(t), p(t))
    if zeta == 0
        % chosen mode has to be 4
        disp('no task request generated!')
        chosen_mode(t) = 4;
    else
        % chosen_mode is chosen from {1, 2, 3}
        disp('task request generated!')
        % task request exists, generate the channel power gain
        h = exprnd(g0 / power(d, 4));
    
        %% step 2.1: solve the optimization problem $\mathcal{P}_{ME}$ (f(t) > 0)
        % calculate f_L and f_U
        f_L = max(sqrt(E_min / (k * W)), W / tau_d);
        f_U = min(sqrt(E_max / (k * W)), f_max);
        if f_L <= f_U
            % the sub-problem is feasible
            disp('mobile execution ($\mathcal{P}_{ME}$) is feasible!')
            
            if B_hat(t) < 0
                f_0 = (V / (-2 * B_hat(t) * k))^(1/3);
            else
                % complex number may exist, which may lead to error
                f_0 = -(V / (2 * B_hat(t) * k))^(1/3);
            end
            
            if (f_0 > f_U && B_hat(t) < 0) || (B_hat(t) >= 0)
                f(t) = f_U;
            elseif f_0 >= f_L && f_0 <= f_U && B_hat(t) < 0
                f(t) = f_0;
            elseif f_0 < f_L && B_hat(t) < 0
                f(t) = f_L;
            end
            % check whether f(t) is zero
            if f(t) == 0
                disp('Something wrong! f is 0!')
            end
            
            % calculate the delay of mobile execution
            cost(t, 1) = W / f(t);
            % calculate the energy consumption of mobile execution
            E(t, 1) = k * W * (f(t)^2);
            % calculate the value of optimization goal
            J_m = -B_hat(t) * k * W * (f(t))^2 + V * W / f(t);
        else
            % the sub-problem is not fasible because (i) the limited 
            % computation capacity or (ii) time cosumed out of deadline or 
            % (iii) the energy consumed out of battery energy level
            % If it is not feasible, it just means that we cannot choose 
            % 'I_m=1'. It dosen't mean that the task has to be dropped.
            disp('mobile execution ($\mathcal{P}_{ME}$) is not feasible!')
            f(t) = 0;
            cost(t, 1) = 0;
            E(t, 1) = 0;
            % 'I_m=1' can never be chosen if mobile execution goal is inf
            J_m = inf;
        % Attention! We do not check whether the energy cunsumed is larger than
        % battery energy level because the problem $\mathcal{J}_{CO}$ does
        % not have constraint (8).
        end
        
        %% step 2.2: solve the optimization problem $\mathcal{P}_{SE}$ (p(t) > 0)
        E_tmp = sigma * L * log(2) / (omega * h);
        p_L_taud = (power(2, L / (omega * tau_d)) - 1) * sigma / h;
        % calculate p_L
        if E_tmp >= E_min
            p_L = p_L_taud;
        else
            % calculate p_E_min (use inline function and fsolve)
            y = @(x) x * L - omega * log2(1 + h*x/sigma) * E_min;
            % accroding to the function figure, p_L_taud is a positive 
            % number around 0.2
            p_E_min = fsolve(y, 0.2, opt);
            p_L = max(p_L_taud, p_E_min);
        end
        % calculate p_U
        if E_tmp >= E_max
            p_U = 0;
        else
            % caculate p_E_max (use inline function and fsolve)
            y = @(x) x * L - omega * log2(1 + h*x/sigma) * E_max;
            % accroding to the function figure, p_E_max is a large positive
            % number around 20
            p_E_max = fsolve(y, 100, opt);
            p_U = min(p_tx_max, p_E_max);
        end
        
        if p_L <= p_U
            % the sub-problem is feasible
            disp('MEC server execution ($\mathcal{P}_{SE}$) is feasible!')
            % calculate p_0
            virtual_battery = B_hat(t);
            y = @(x) virtual_battery * log2(1 + h*x/sigma) + ...
                h * (V - virtual_battery*x) / log(2) / (sigma + h*x);
            p_0 = fsolve(y, 0.5, opt);

            if (p_U < p_0 && B_hat(t) < 0) || B_hat(t) >= 0
                p(t) = p_U;
            elseif p_0 < p_L && B_hat(t) < 0
                p(t) = p_L;
            elseif p_0 >= p_L && p_0 <= p_U && B_hat(t) < 0
                p(t) = p_0;
            end
            % check whether p(t) is zero
            if p(t) == 0
                disp('Something wrong! p is 0!')
            end
            
            % calculate the delay of MEC server execution
            cost(t, 2) = L / (omega * log2(1 + h*p(t)/sigma));
            % calculate the energy consumption of MEC server execution
            E(t, 2) = p(t) * cost(t, 2);
            % calculate the value of optimization goal
            J_s = (-B_hat(t) * p(t) + V) * cost(t, 2);
        else
            % the sub-problem is not feasible because (i) the limited transmit 
            % power or (ii) time cosumed out of deadline or (iii) the energy 
            % consumed out of battery energy level
            % If it is not feasible, it just means that we cannot choose 
            % 'I_s=1'. It dosen't mean that the task has to be dropped.
            disp('MEC server execution ($\mathcal{P}_{SE}$) is not feasible!')
            p(t) = 0;
            cost(t, 2) = 0;
            E(t, 2) = 0;
            % 'I_s=1' can never be chosen if MEC server execution goal is inf
            J_s = inf;
        % Similarly, we do not check whether the energy cunsumed is larger than
        % battery energy level because the problem $\mathcal{J}_{CO}$ does
        % not have constraint (8).
        end
        
        %% step 3: choose the best execution mode
        J_d = V * phi;
        disp(['J_m:', num2str(J_m)])
        disp(['J_s:', num2str(J_s)])
        [~, mode] = min([J_m, J_s, J_d]);
        chosen_mode(t) = mode;
    end
    
    %% step 4: according to the chosen execution mode, calculate the real dealy and energy consumption
    if chosen_mode(t) == 1
        % mobile execution is chosen
        cost(t, 3) = cost(t, 1);
        E(t, 3) = E(t, 1);
    elseif chosen_mode(t) == 2
        % MEC server execution is chosen
        cost(t, 3) = cost(t, 2);
        E(t, 3) = E(t, 2);
    elseif chosen_mode(t) == 3
        % task is dropped, the delay is the task dropping penalty and the 
        % energy consumption is zero
        cost(t, 3) = phi;
        E(t, 3) = 0;
    else
        % no task is requested, the delay and the energy consumption are
        % both zero
        cost(t, 3) = 0;
        E(t, 3) = 0;
    end
    
    %% step 5: update the battery energy level and go to next time slot
    B(t + 1) = B(t) - E(t, 3) + e(t);
    t = t + 1;
end

%% step 6: evaluate the simulation results
% 1. the battery energy level vs. time slot
figure
plot(1:T, B(1:T));
hold on
plot(1:T, repmat(theta + E_H_max, [T, 1]), '-')
title('Envolution of battery energy level')
xlabel('time slot')
ylabel('battery energy level $B_t$', 'Interpreter','latex')

% 2. the average execution cost vs. time slot
accumulated = 0;
average_cost = zeros(T, 1);
request_num = 0;
for t = 1: T
    accumulated = accumulated + cost(t, 3);
    if cost(t, 3) ~= 0
        % there exists task request
        request_num = request_num + 1;
    end
    average_cost(t) = accumulated / request_num;
end
figure
plot(1:T, average_cost);
title('Envolution of average execution cost')
xlabel('time slot')
ylabel('average execution cost $\frac{1}{T} \sum_{t=0}^{T-1} cost^t$', 'Interpreter','latex')

% 3. the average ratio of each chosen mode vs. time slot
average_ratio = zeros(T, 3);
mobile_exe = 0; server_exe = 0; drop = 0;
request_num = 0;
for t = 1: T
    if cost(t, 3) == 0
        continue
    else
        request_num = request_num + 1;
        if chosen_mode(t) == 1
            mobile_exe = mobile_exe + 1;
        elseif chosen_mode(t) == 2
            server_exe = server_exe + 1;
        else
            drop = drop + 1;
        end
    end
    average_ratio(t, :) = [mobile_exe, server_exe, drop] / request_num;
end
figure
plot(1:T, average_ratio(:, 1));
hold on
plot(1:T, average_ratio(:, 2));
hold on
plot(1:T, average_ratio(:, 3));
legend('mobile execution', 'MEC server execution', 'drop')
title('Envolution of average ratio of chosen modes')
xlabel('time slot')
ylabel('average  ratio of chosen modes $\frac{1}{T} \sum_{t=0}^{T-1} \{I_m^t, I_s^t, I_d^t\}$', 'Interpreter','latex')

%% This script simulates the LODCO-based Greedy algorithm.
%  author: Hailiang Zhao
clc, clear
opt = optimset('Display', 'none');

%% basic parameter settings
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
d0 = 1;                   % the relative distance between each mobile device and each MEC server

%% parameter control
N = 20;                   % the number of mobile devices
M = 8;                    % the number of MEC servers
T = 500;                  % the number of time slot (a.k.a. the size of the time horizon)
tau_d = 0.002;            % execution deadline (in second)
d = 50;                   % the distance between the mobile device and the MEC server (in meter)
E_min = 0.02e-3;          % the minimum amout of battery output energy (in J)
V = 1e-5;                 % the weight of penalty (the control parameter introduced by Lyapunov Optimization)
rho = 0.6;                % the probability that the computation task is requested
max_connects = 4;         % the maximum number of processible mobile devices for each MEC server ($ \frac{f_s^{max} \tau}{L X} $)

% the lower bound of perturbation parameter
E_max_hat = min(max(k * W * (f_max)^2, p_tx_max * tau), E_max);
theta = E_max_hat + V * phi / E_min;

%% allocate storage for valuable results
B = zeros(T, N);                   % the battery energy level (in J)
B_hat = zeros(T, N);               % the virtual battery energy level ($B_hat = B - theta$)
e = zeros(T, N);                   % the amout of the harvested and stored energy (in J)
chosen_mode = zeros(T, N);         % {1: local, 2: remote, 3: drop, 4: no task request}
chosen_server = zeros(T, N);       % record the index of chosen server for each mobile device if its choice is MEC server execution
f = zeros(T, N);                   % the CPU-cycle frequency of local execution (in Hz)
p = zeros(T, N);                   % the transmit power of computation offloading (in W)

mobile_exe_cost = zeros(T, N);     % the mobile execution cost (delay) (in second)
server_exe_cost = zeros(T, N);     % the MEC server execution cost (delay) (in second)
final_chosen_cost = zeros(T, N);   % the final execution delay under currently chosen modes (in second)

mobile_exe_E = zeros(T, N);        % the energy consumption for mobile execution (in J)
server_exe_E = zeros(T, N);        % the energy consumption for MEC server execution (in J)
final_chosen_E = zeros(T, N);      % the energy consumption of the final chosen modes (in J)

%% simulation begin
t = 1;
while t <= T
    disp(['===> Time slot #', num2str(t), ' <==='])
    
    %% allocate storage for mode-chosen
    device_server_pairs = [];      % each column represents i, j, J_s^{\star}(i, j), respectively
    remained_connects = max_connects * ones(M, 1);    % the available connections of MEC servers
    J_m = zeros(N, 1); J_s = zeros(N, M);             % the matrices for J_m and J_s values
    p_mat = zeros(N, M);                              % the matrix for transmit power
    server_cost_mat = zeros(N, M);                    % the matrix for MEC server execution cost
    server_E_mat = zeros(N, M);                       % the matrix for energy consumption of MEC server
    
    %% initialization
    % generate the virtual battery energy level
    B_hat(t, :) = B(t, :) - theta;
    % generate the channel power gain (from each mobile device to each MEC sever)
    distances = unifrnd(20, 80, N, M);
    gamma = exprnd(1, N, M);
    h_mat = g0 * gamma .* power(d0 ./ distances, 4);
    
    %% step 1: for each mobile device, choose the initial mode
    for i = 1: N
        disp(['Mobile device #', num2str(i)])
        
        %% step 1.1: get the optimal energy harvesting no matter whether task is requested
        E_H_t = unifrnd(0, E_H_max);
        if B_hat(t, i) <= 0
            e(t, i) = E_H_t;
        end
        
        %% step 1.2: get the (initial) optimal computation offloading strategy (I_m(i), I_s(i, :), I_d(i), f(t, i), p(t, i))
        % generate the task request
        zeta = binornd(1, rho);
        if zeta == 0
            % chosen mode has to be 4
            disp('no task request generated!')
            chosen_mode(t, i) = 4;
        else
            % chosen_mode is chosen from {1, 2, 3}
            disp('task request generated!')
            
            %% step 1.2.1: solve the optimization problem $\mathcal{P}_{ME}$ (f(t, i) > 0)
            % calculate f_L and f_U
            f_L = max(sqrt(E_min / (k * W)), W / tau_d);
            f_U = min(sqrt(E_max / (k * W)), f_max);
            if f_L <= f_U
                % the sub-problem is feasible
                disp('mobile execution ($\mathcal{P}_{ME}$) is feasible!')
                
                if B_hat(t, i) < 0
                    f_0 = (V / (-2 * B_hat(t, i) * k))^(1/3);
                else
                    % complex number may exist, which may lead to error
                    f_0 = -(V / (2 * B_hat(t, i) * k))^(1/3);
                end
                
                if (f_0 > f_U && B_hat(t, i) < 0) || (B_hat(t, i) >= 0)
                    f(t, i) = f_U;
                elseif f_0 >= f_L && f_0 <= f_U && B_hat(t, i) < 0
                    f(t, i) = f_0;
                elseif f_0 < f_L && B_hat(t, i) < 0
                    f(t, i) = f_L;
                end
                % check whether f(t, i) is zero
                if f(t, i) == 0
                    disp('Something wrong! f is 0!')
                end
                
                % calculate the delay of mobile execution
                mobile_exe_cost(t, i) = W / f(t, i);
                % calculate the energy consumption of mobile execution
                mobile_exe_E(t, i) = k * W * (f(t, i)^2);
                % calculate the value of optimization goal
                J_m(i) = -B_hat(t, i) * k * W * (f(t, i)^2 + V * W / f(t, i));
            else
                % the sub-problem is not fasible because (i) the limited 
                % computation capacity or (ii) time cosumed out of deadline 
                % or (iii) the energy consumed out of battery energy level
                % If it is not feasible, it just means that we cannot choose 
                % 'I_m(i)=1'. It dosen't mean that the task has to be dropped.
                disp('mobile execution ($\mathcal{P}_{ME}$) is not feasible!')
                f(t, i) = 0;
                mobile_exe_cost(t, i) = 0;
                mobile_exe_E(t, i) = 0;
                % 'I_m(i)=1' can never be chosen if mobile execution goal is inf
                J_m(i) = inf;
            end
            
            %% step 1.2.2: solve the optimization problem $\mathcal{P}_{SE}$ (p(t, i) > 0)
            % calculate J_s(i, j) from mobile device i to each MEC server j
            for j = 1: M
                disp(['MEC server #', num2str(j)])
                h = h_mat(i, j);
                
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
                    virtual_battery = B_hat(t, i);
                    y = @(x) virtual_battery * log2(1 + h*x/sigma) + ...
                        h * (V - virtual_battery*x) / log(2) / (sigma + h*x);
                    p_0 = fsolve(y, 0.5, opt);
                    
                    if (p_U < p_0 && B_hat(t, i) < 0) || B_hat(t, i) >= 0
                        p_mat(i, j) = p_U;
                    elseif p_0 < p_L && B_hat(t) < 0
                        p_mat(i, j) = p_L;
                    elseif p_0 >= p_L && p_0 <= p_U && B_hat(t) < 0
                        p_mat(i, j) = p_0;
                    end
                    % check whether p_mat(i, j) is zero
                    if p_mat(i, j) == 0
                        disp('Something wrong! p is 0!')
                    end
                    
                    % calculate the delay of MEC server execution
                    server_cost_mat(i, j) = L / (omega * log2(1 + h*p_mat(i, j)/sigma));
                    % calculate the energy consumption of MEC server execution
                    server_E_mat(i, j) = p_mat(i, j) * server_cost_mat(i, j);
                    % calculate the value of optimization goal
                    J_s(i, j) = (-B_hat(t, i) * p_mat(i, j) + V) * server_cost_mat(i, j);
                    
                    % (we can not set server_exe_cost(t, i) and server_exe_E(t, i) for now)
                else
                    % the sub-problem is not feasible because (i) the limited transmit 
                    % power or (ii) time cosumed out of deadline or (iii) the energy 
                    % consumed out of battery energy level
                    % If it is not feasible, it just means that we cannot choose 
                    % 'I_s(i,j)=1'. It dosen't mean that the task has to be dropped.
                    disp('MEC server execution ($\mathcal{P}_{SE}$) is not feasible!')
                    p_mat(i, j) = 0;
                    server_cost_mat(i, j) = 0;
                    server_E_mat(i, j) = 0;
                    % 'I_s(i,j)=1' can never be chosen if MEC server execution goal is inf
                    J_s(i, j) = inf;
                    
                    % (we can not set server_exe_cost(t, i) and server_exe_E(t, i) for now)
                % Similarly, we do not check whether the energy cunsumed is larger than
                % battery energy level because the problem $\mathcal{J}_{CO}$ does
                % not have constraint (8).
                end
            end
            
            %% step 1.2.3: choose the (initial) optimal execution mode
            J_d = V * phi;
            disp(['J_m(i):', num2str(J_m(i))])
            disp(['J_s(i,:)', num2str(J_s(i, :))])
            [~, mode] = min([J_m(i), J_s(i, :), J_d]);
            if mode == 1
                % mobile execution
                chosen_mode(t, i) = 1;
                final_chosen_cost(t, i) = mobile_exe_cost(t, i);
                final_chosen_E(t, i) = mobile_exe_E(t, i);
            elseif mode == (M+2)
                % drop
                chosen_mode(t, i) = 3;
                final_chosen_cost(t, i) = phi;
                final_chosen_E(t, i) = 0;
            else
                % MEC servre execution
                chosen_mode(t, i) = 2;
                % add i, the chosen j, and their J_s(i, j) value into the pairs
                device_server_pairs = [device_server_pairs; [i, mode-1, J_s(i, mode-1)]];
            end
        end
    end
    
    %% step 2: allocate connection right for each mobile device who chooses MEC server execution
    while ~isempty(device_server_pairs)
        for j = 1: M
            % find every pair who chooses the MEC server j
            device_j_pairs = device_server_pairs(device_server_pairs(:, 2) == j, :);
            % find those devices is
            is = device_j_pairs(:, 1);
            if isempty(is)
                disp(['For current MEC server #', num2str(j), ', no mobile device choose it!'])
                % go to handle next MEC server
                continue;
            end
            if remained_connects(j) >= length(is)
                % every mobile device who chooses j can be connected,
                % set their final modes as 2 and record the final cost and energy consumption for them
                chosen_mode(t, is) = 2;     % this is not necessary
                p(t, is) = transp(p_mat(is, j));                            % those assignment might exist problem
                server_exe_cost(t, is) = transp(server_cost_mat(is, j));    % those assignment might exist problem
                server_exe_E(t, is) = transp(server_E_mat(is, j));          % those assignment might exist problem
                chosen_server(t, is) = repmat(j, 1, length(is));            % those assignment might exist problem
                final_chosen_cost(t, is) = server_exe_cost(t, is);          % those assignment might exist problem
                final_chosen_E(t, is) = server_exe_E(t, is);                % those assignment might exist problem
                
                % update remained_connects for j
                remained_connects(j) = remained_connects(j) - length(is);
                % finally, remove them from global pairs, they are done
                device_server_pairs(device_server_pairs(:, 2) == j, :) = [];
            else
                if remained_connects(j) == 0
                    % no mobile device can be connected, remove all mobile devices who chooses j from global pairs
                    device_server_pairs(device_server_pairs(:, 2) == j, :) = [];
                    % set those mobile devices' J_s(i, j) as inf
                    % (which means no matter what happens, the MEC server j can never be chosen)
                    J_s(is, j) = inf;
                    % choose mode again for those i and insert new potential pairs into global pairs
                    for idx = 1: numel(is)
                        i = is(idx);
                        [~, mode] = min([J_m(i), J_s(i, :), J_d]);
                        if mode == 1
                            chosen_mode(t, i) = 1;
                        elseif mode == (M+2)
                            chosen_mode(t, i) = 3;
                        else
                            chosen_mode(t, i) = 2;
                            device_server_pairs = [device_server_pairs; [i, mode-1, J_s(i, mode-1)]];
                        end
                    end
                else
                    % some mobile devices can be connected, set their final modes as 2 and remove them from global pairs
                    % besides, record the final cost and energy consumption for them
                    % for the left i', remove them from global pairs and set J_s(i', j) as inf, 
                    % choose mode again for those i' and insert new potential pairs into global pairs
                    
                    % sort the J_s(is, j) and return those lucky idxs
                    [~, idxs] = sort(device_j_pairs(:, 3));
                    for idx = 1: remained_connects(j)
                        i = idxs(idx);
                        chosen_mode(t, i) = 2;
                        p(t, i) = p_mat(i, j);
                        server_exe_cost(t, i) = server_cost_mat(i, j);
                        server_exe_E(t, i) = server_E_mat(i, j);
                        chosen_server(t, i) = j;
                        final_chosen_cost(t, i) = server_exe_cost(t, i);
                        final_chosen_E(t, i) = server_exe_E(t, i);
                        
                        % remove i from global pairs
                        device_server_pairs(device_server_pairs(:, 1) == i, :) = [];
                    end
                    
                    % update remained_connects for j
                    remained_connects(j) = 0;
                    
                    % for those unlucky mobile devices i', set J_s(i', j) as inf (i' are in currently device_server_pairs(: , j) now)
                    residual_is = device_server_pairs(device_server_pairs(:, 2) == j, 1);
                    J_s(residual_is, j) = inf;
                    for idx = 1: numel(residual_is)
                        residual_i = residual_is(idx);
                        [~, mode] = min([J_m(residual_i), J_s(residual_i, :), J_d]);
                        if mode == 1
                            chosen_mode(t, residual_i) = 1;
                        elseif mode == (M+2)
                            chosen_mode(t, residual_i) = 3;
                        else
                            chosen_mode(t, residual_i) = 2;
                            device_server_pairs = [device_server_pairs; [residual_i, mode-1, J_s(residual_i, mode-1)]];
                        end
                    end
                end
            end
        end
    end
    
    %% step 3: update the battery energy level and go to the next time slot
    B(t + 1, :) = B(t, :) - final_chosen_E(t, :) + e(t, :);
    t = t + 1;
    
end

%% step 4: evaluate the simulation results
% 1. the battery energy level vs. time slot
figure
plot(1:T, B(1:T, :))
hold on
plot(1:T, repmat(theta + E_H_max, [T, 1]), '-')
title('Envolution of battery energy level')
xlabel('time slot')
ylabel('battery energy level $B_t$ of each mobile device', 'Interpreter','latex')

% 2. the average execution cost vs. time slot
accumulated = 0;
average_cost = zeros(T, N);
figure
for i = 1: N
    % draw for each mobile device
    request_num = 0;
    for t = 1: T
        accumulated = accumulated + final_chosen_cost(t, i);
        if chosen_mode(t, i) ~= 4
            % there exists task request
            request_num = request_num + 1;
        end
        average_cost(t, i) = accumulated / request_num;
    end
    plot(1:T, average_cost(:, i));
end
title('Envolution of average execution cost')
xlabel('time slot')
ylabel('average execution cost $\frac{1}{T} \sum_{t=0}^{T-1} cost^t$ of each mobile device', 'Interpreter','latex')

% 3. the average ratio of each chosen mode of the ith mobile device vs. time slot
average_ratio = zeros(T, 3);
mobile_exe = 0; server_exe = 0; drop = 0;
request_num = 0;
i = 1;          % we simply choose the first mobile device
for t = 1: T
    if final_chosen_cost(t, i) == 0
        continue
    else
        request_num = request_num + 1;
        if chosen_mode(t, i) == 1
            mobile_exe = mobile_exe + 1;
        elseif chosen_mode(t, i) == 2
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
ylabel('average  ratio of chosen modes $\frac{1}{T} \sum_{t=0}^{T-1} \{I_m^t, I_s^t, I_d^t\}$ of the i-th mobile device', 'Interpreter','latex')


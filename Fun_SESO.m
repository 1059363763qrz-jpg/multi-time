function [P_seso_to_dso_charge, P_seso_to_dso_discharge, P_seso_from_mg1_c, P_seso_from_mg1_d, P_seso_from_mg2_c, P_seso_from_mg2_d, P_seso_from_mg3_c, P_seso_from_mg3_d, obj_seso, E_seso_val, P_seso_ch_val, P_seso_dis_val] = Fun_SESO(...
    P_dso_charge, P_dso_discharge, lambda_dso_charge, lambda_dso_discharge, ...
    P_mg1_lease_c, P_mg1_lease_d, lambda_mg1_c, lambda_mg1_d, ...
    P_mg2_lease_c, P_mg2_lease_d, lambda_mg2_c, lambda_mg2_d, ...
    P_mg3_lease_c, P_mg3_lease_d, lambda_mg3_c, lambda_mg3_d, rho)
%Fun_SESO: 共享储能运营商模型

%% 1. 参数
T = 24;
% SESO自有储能参数
SOC_min = 0.1;
SOC_max = 0.9;
E_max = 2000; % 自有储能容量(kWh)，设置为单独无法满足配网需求
P_ch_max = 500; % 最大充电功率(kW)
P_dis_max = 500; % 最大放电功率(kW)
eta_ch = 0.95; eta_dis = 0.95;
initial_SOC = 0.5;

% 交互价格
price_dso_charge = 0.25;  % 向DSO出租充电服务的价格
price_dso_discharge = 0.55; % 向DSO购买放电服务的价格
price_mg_lease = 0.22; % 向MG租赁储能服务的价格

%% 2. 决策变量
% 自有储能变量
E_seso = sdpvar(1, T+1); P_seso_ch = sdpvar(1, T); P_seso_dis = sdpvar(1, T);
u_seso_ch = binvar(1,T); u_seso_dis = binvar(1,T);

% 从微网租赁的功率
P_seso_from_mg1_c = sdpvar(1, T); P_seso_from_mg1_d = sdpvar(1, T);
P_seso_from_mg2_c = sdpvar(1, T); P_seso_from_mg2_d = sdpvar(1, T);
P_seso_from_mg3_c = sdpvar(1, T); P_seso_from_mg3_d = sdpvar(1, T);

% 提供给DSO的功率
P_seso_to_dso_charge = sdpvar(1, T);
P_seso_to_dso_discharge = sdpvar(1, T);

%% 3. 约束条件
C = [];

% 3.1 SESO自有储能约束
C = [C, E_seso(1) == initial_SOC * E_max];
for t = 1:T
    C = [C, E_seso(t+1) == E_seso(t) + P_seso_ch(t) * eta_ch - P_seso_dis(t) / eta_dis];
    C = [C, SOC_min * E_max <= E_seso(t+1) <= SOC_max * E_max];
    C = [C, 0 <= P_seso_ch(t) <= u_seso_ch(t) * P_ch_max];
    C = [C, 0 <= P_seso_dis(t) <= u_seso_dis(t) * P_dis_max];
    C = [C, u_seso_ch(t) + u_seso_dis(t) <= 1];
end
C = [C, E_seso(T) == initial_SOC * E_max];

% 3.2 SESO功率池平衡约束
for t = 1:T
    % 充电功率池：从各微网租来的充电服务 + 自有储能充电 = 提供给DSO的总充电服务
    total_charge_in = P_seso_from_mg1_c(t) + P_seso_from_mg2_c(t) + P_seso_from_mg3_c(t) + P_seso_ch(t);
    C = [C, total_charge_in == P_seso_to_dso_charge(t)];
    
    % 放电功率池：从各微网租来的放电服务 + 自有储能放电 = 提供给DSO的总放电服务
    total_discharge_in = P_seso_from_mg1_d(t) + P_seso_from_mg2_d(t) + P_seso_from_mg3_d(t) + P_seso_dis(t);
    C = [C, total_discharge_in == P_seso_to_dso_discharge(t)];
end
%3-3交互变量约束
for t=1:T
   C=[C,0<=P_seso_from_mg1_c<=500,];
   C=[C,0<=P_seso_from_mg2_c<=500,];
   C=[C,0<=P_seso_from_mg3_c<=500,];
   C=[C, 0 <= P_seso_from_mg1_d(t) <= 500]; 
   C=[C, 0 <= P_seso_from_mg2_d(t) <= 500];
   C=[C, 0 <= P_seso_from_mg3_d(t) <= 500];
   C=[C, 0 <= P_seso_to_dso_charge(t) <= 1000];
   C=[C, 0 <= P_seso_to_dso_discharge(t) <= 1000];
end
%% 4. 目标函数 (最大化利润)
% 来自DSO的收益
revenue_dso = sum(price_dso_charge .* P_seso_to_dso_charge) - sum(price_dso_discharge .* P_seso_to_dso_discharge);
% 支付给微网的租赁成本+自有储能的折现成本（折现系数设置小于微网租赁系数，即优先使用自有储能）
cost = sum(price_mg_lease * (P_seso_from_mg1_c + P_seso_from_mg1_d + P_seso_from_mg2_c + P_seso_from_mg2_d + P_seso_from_mg3_c + P_seso_from_mg3_d))+sum(0.01*(P_seso_ch+P_seso_dis));

% ADMM罚项
admm_penalty = (rho/2) * (norm(P_seso_to_dso_charge - P_dso_charge)^2 + norm(P_seso_to_dso_discharge - P_dso_discharge)^2 ...
            + norm(P_seso_from_mg1_c - P_mg1_lease_c)^2 + norm(P_seso_from_mg1_d - P_mg1_lease_d)^2 ...
            + norm(P_seso_from_mg2_c - P_mg2_lease_c)^2 + norm(P_seso_from_mg2_d - P_mg2_lease_d)^2 ...
            + norm(P_seso_from_mg3_c - P_mg3_lease_c)^2 + norm(P_seso_from_mg3_d - P_mg3_lease_d)^2) ...
            - sum(lambda_dso_charge .* (P_seso_to_dso_charge - P_dso_charge)) ...
            - sum(lambda_dso_discharge .* (P_seso_to_dso_discharge - P_dso_discharge)) ...
            - sum(lambda_mg1_c .* (P_seso_from_mg1_c - P_mg1_lease_c)) - sum(lambda_mg1_d .* (P_seso_from_mg1_d - P_mg1_lease_d)) ...
            - sum(lambda_mg2_c .* (P_seso_from_mg2_c - P_mg2_lease_c)) - sum(lambda_mg2_d .* (P_seso_from_mg2_d - P_mg2_lease_d)) ...
            - sum(lambda_mg3_c .* (P_seso_from_mg3_c - P_mg3_lease_c)) - sum(lambda_mg3_d .* (P_seso_from_mg3_d - P_mg3_lease_d));
            
% 目标：最小化 (负的利润)
Objective = -(revenue_dso - cost) + admm_penalty;

%% 5. 求解与输出
ops = sdpsettings('solver','gurobi','verbose',0','debug', 0);
ops.cplex.parameters.output.clock = 0;  % 关闭时钟输出
ops.cplex.parameters.output.solution = 0;  % 关闭解信息输出
ops.cplex.display = 'off';
solution = solvesdp(C, Objective, ops);
if solution.problem ~= 0
    % 如果problem不为0，说明求解失败
    fprintf('!!! 警告: 在 Fun_seso 中求解失败 !!!\n');
    disp(yalmiperror(solution.problem)); % 打印具体的错误信息
end
P_seso_to_dso_charge = value(P_seso_to_dso_charge);
P_seso_to_dso_discharge = value(P_seso_to_dso_discharge);
P_seso_from_mg1_c = value(P_seso_from_mg1_c); P_seso_from_mg1_d = value(P_seso_from_mg1_d);
P_seso_from_mg2_c = value(P_seso_from_mg2_c); P_seso_from_mg2_d = value(P_seso_from_mg2_d);
P_seso_from_mg3_c = value(P_seso_from_mg3_c); P_seso_from_mg3_d = value(P_seso_from_mg3_d);
obj_seso = value(sum(0.01*(P_seso_ch+P_seso_dis)));
E_seso_val = value(E_seso);
P_seso_ch_val = value(P_seso_ch);
P_seso_dis_val = value(P_seso_dis);
end
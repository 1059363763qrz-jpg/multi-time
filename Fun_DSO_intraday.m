function [delta_charge, delta_discharge, delta_grid] = Fun_DSO_intraday(plan_DA, dso_data_rt, k, H, ...
    delta_seso_c_k, delta_seso_d_k, lambda_c, lambda_d, delta_mg1_net, delta_mg2_net, delta_mg3_net, rho)

    yalmip('clear');
    T = H; % 优化时域为滚动窗口长度

    %% 权重系数
    w_adjust_grid = 0.005;    % 主网功率调整成本系数
    w_adjust_dr = 0.002;     % 需求响应调整成本
    w_adjust_ses = 0.003;
    w_unbalance = 100;    % 功率不平衡惩罚系数

    %% 决策变量 (调整量)
    delta_charge = sdpvar(1, T);
    delta_discharge = sdpvar(1, T);
    delta_grid = sdpvar(1, T);
    delta_Lshift = sdpvar(1, T);
    delta_SIL = sdpvar(1, T);
    P_unbalance = sdpvar(1, T); % 松弛变量

    %% 目标函数
    % 1. Minimize power adjustment costs
     cost_adjust = w_adjust_grid*sum(delta_grid.^2)+w_adjust_dr*sum(delta_Lshift.^2+delta_SIL.^2)+w_adjust_ses*sum(delta_charge.^2+delta_discharge.^2);
    
    % 2. Minimize power unbalance penalty (Linearized)
    % MODIFICATION: Changed from sum(P_unbalance.^2) to sum(abs(P_unbalance))
    cost_unbalance = w_unbalance * sum(P_unbalance.^2);

    % 3. ADMM penalty term
    admm_penalty = sum(lambda_c .* (delta_charge - delta_seso_c_k)) ...
                 + sum(lambda_d .* (delta_discharge - delta_seso_d_k)) ...
                 + (rho/2) * (norm(delta_charge - delta_seso_c_k)^2 + norm(delta_discharge - delta_seso_d_k)^2);

     Objective = cost_adjust+cost_unbalance+admm_penalty;
    % Objective =  admm_penalty;
      % Objective = cost_adjust  + admm_penalty;
      % Objective =  cost_unbalance + admm_penalty; 

    %% 约束条件
    C = [];
    for t = 1:T
        % 获取对应的日前小时索引
        hour_idx = floor(((k+t-1)-1) * 15 / 60) + 1;
        
        % 实际功率 = 日前计划 + 调整量
        P_charge_actual = plan_DA.P_dso_charge(hour_idx) + delta_charge(t);
        P_discharge_actual = plan_DA.P_dso_discharge(hour_idx) + delta_discharge(t);
        P_grid_actual = plan_DA.P_dso_grid(hour_idx) + delta_grid(t);
        P_mg1_net_actual = plan_DA.P_mg1_net_exchange(hour_idx) + delta_mg1_net(t);
        P_mg2_net_actual = plan_DA.P_mg2_net_exchange(hour_idx) + delta_mg2_net(t);
        P_mg3_net_actual = plan_DA.P_mg3_net_exchange(hour_idx) + delta_mg3_net(t);
        Lshift_actual = plan_DA.Lshift_effect_dso(hour_idx) + delta_Lshift(t);
        SIL_actual = plan_DA.SIL_effect_dso(hour_idx) + delta_SIL(t);
        % ... (MG2, MG3)
        
        % 功率平衡约束 (简化模型)
        % 总注入 = 总流出
        total_gen = P_grid_actual + sum(dso_data_rt.P_pv(:,t)) + sum(dso_data_rt.P_wt(:,t)) + P_discharge_actual+ P_mg1_net_actual...
                     +P_mg2_net_actual+P_mg3_net_actual;
        total_load = sum(dso_data_rt.PLOAD(:,t)) + Lshift_actual - SIL_actual+P_charge_actual ; %... + MG2, MG3

        % 此处网损用日前计划值近似
        P_loss_planned = plan_DA.P_loss(hour_idx); % 简化为0，或从plan_DA中提取
        
        C = [C, total_gen == total_load + P_loss_planned + P_unbalance(t)];
        
        % 调整量边界约束 (防止调整过大)
        C = [C, -500 <= delta_grid <= 500]; % 示例
        C = [C, P_charge_actual >= 0, P_discharge_actual >= 0];
        C = [C, P_charge_actual <=1000, P_discharge_actual <= 1000];
        C = [C, -0.05*sum(dso_data_rt.PLOAD(:,t)) <= delta_Lshift(t) <= 0.05*sum(dso_data_rt.PLOAD(:,t))];
        C = [C, -0.05*sum(dso_data_rt.PLOAD(:,t)) <= delta_SIL(t) <= 0.05*sum(dso_data_rt.PLOAD(:,t))];
    end
    
    %% 求解
    ops = sdpsettings('solver','cplex','verbose',0);
    sol = optimize(C, Objective, ops);
    
    if sol.problem == 0
        delta_charge = value(delta_charge);
        delta_discharge = value(delta_discharge);
        delta_grid = value(delta_grid);
    else
        disp('!!! DSO 日内子问题求解失败 !!!');
        delta_charge = zeros(1,T); delta_discharge = zeros(1,T); delta_grid = zeros(1,T);
    end
end
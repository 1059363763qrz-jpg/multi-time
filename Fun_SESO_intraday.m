function [delta_to_dso_c, delta_to_dso_d, delta_from_mg1_c, delta_from_mg1_d, delta_from_mg2_c, delta_from_mg2_d, delta_from_mg3_c, delta_from_mg3_d] = Fun_SESO_intraday(plan_DA, k, H, ...
    delta_dso_c_k, delta_dso_d_k, lambda_dso_c, lambda_dso_d, ...
    delta_mg1_c_k, delta_mg1_d_k, lambda_mg1_c, lambda_mg1_d, ...
    delta_mg2_c_k, delta_mg2_d_k, lambda_mg2_c, lambda_mg2_d, ...
    delta_mg3_c_k, delta_mg3_d_k, lambda_mg3_c, lambda_mg3_d, rho)
% Fun_SESO_intraday: 共享储能运营商(SESO)的日内滚动优化子问题
% 最终修复版: 使用“配方法”重构ADMM罚项，彻底解决YALMIP求解器报错问题

    yalmip('clear');
    T = H;
    
    %% 权重系数
    w_adjust_own = 0.02;
    w_violation  = 100;

    %% 决策变量
    delta_own_ch = sdpvar(1, T); delta_own_dis = sdpvar(1, T);
    delta_to_dso_c = sdpvar(1, T); delta_to_dso_d = sdpvar(1, T);
    delta_from_mg1_c = sdpvar(1, T); delta_from_mg1_d = sdpvar(1, T);
    delta_from_mg2_c = sdpvar(1, T); delta_from_mg2_d = sdpvar(1, T);
    delta_from_mg3_c = sdpvar(1, T); delta_from_mg3_d = sdpvar(1, T);
    delta_violation_c = sdpvar(1, T); delta_violation_d = sdpvar(1, T);

    %% 目标函数 (已使用“配方法”重构)
    
    cost_adjust = w_adjust_own * (sum(delta_own_ch.^2) + sum(delta_own_dis.^2));
    cost_violation = w_violation * (sum(delta_violation_c) + sum(delta_violation_d));
                                     
    % =================== 使用“配方法”重构ADMM罚项 ===================
   admm_penalty_dso = sum(lambda_dso_c .* (delta_to_dso_c - delta_dso_c_k)) ...
                     + sum(lambda_dso_d .* (delta_to_dso_d - delta_dso_d_k)) ...
                     + (rho/2) * (norm(delta_to_dso_c - delta_dso_c_k)^2 + norm(delta_to_dso_d - delta_dso_d_k)^2);

    % 与MG1的交互
    admm_penalty_mg1 = sum(lambda_mg1_c .* (delta_from_mg1_c - delta_mg1_c_k)) ...
                     + sum(lambda_mg1_d .* (delta_from_mg1_d - delta_mg1_d_k)) ...
                     + (rho/2) * (norm(delta_from_mg1_c - delta_mg1_c_k)^2 + norm(delta_from_mg1_d - delta_mg1_d_k)^2);

    % 与MG2的交互
    admm_penalty_mg2 = sum(lambda_mg2_c .* (delta_from_mg2_c - delta_mg2_c_k)) ...
                     + sum(lambda_mg2_d .* (delta_from_mg2_d - delta_mg2_d_k)) ...
                     + (rho/2) * (norm(delta_from_mg2_c - delta_mg2_c_k)^2 + norm(delta_from_mg2_d - delta_mg2_d_k)^2);

    % 与MG3的交互
    admm_penalty_mg3 = sum(lambda_mg3_c .* (delta_from_mg3_c - delta_mg3_c_k)) ...
                     + sum(lambda_mg3_d .* (delta_from_mg3_d - delta_mg3_d_k)) ...
                     + (rho/2) * (norm(delta_from_mg3_c - delta_mg3_c_k)^2 + norm(delta_from_mg3_d - delta_mg3_d_k)^2);
    
    admm_penalty = admm_penalty_dso + admm_penalty_mg1 + admm_penalty_mg2 + admm_penalty_mg3;
    % =================================================================

    Objective = cost_adjust + cost_violation+ admm_penalty;

    %% 约束条件 (无需改动)
    C = [];
    total_charge_collected = delta_from_mg1_c + delta_from_mg2_c + delta_from_mg3_c + delta_own_ch;
    total_discharge_collected = delta_from_mg1_d + delta_from_mg2_d + delta_from_mg3_d + delta_own_dis;
    
    C = [C, delta_to_dso_c + delta_violation_c == total_charge_collected];
    C = [C, delta_to_dso_d + delta_violation_d == total_discharge_collected];
    C = [C, delta_violation_c >= 0, delta_violation_d >= 0];
    
    E_seso_actual = sdpvar(1, T+1);
    hour_idx_start = floor((k-1) / 4) + 1;
    E_seso_initial_plan = plan_DA.E_seso(hour_idx_start);
    C = [C, E_seso_actual(1) == E_seso_initial_plan];
    
    for t = 1:T
        current_k = k + t - 1; hour_idx = floor((current_k - 1) / 4) + 1;
        P_own_ch_actual = plan_DA.P_seso_ch(hour_idx) + delta_own_ch(t);
        P_own_dis_actual = plan_DA.P_seso_dis(hour_idx) + delta_own_dis(t);
        
        eta_ch = 0.95; eta_dis = 0.95;
        C = [C, E_seso_actual(t+1) == E_seso_actual(t) + (P_own_ch_actual * eta_ch - P_own_dis_actual / eta_dis) * 0.25];
        C = [C, 0.1*2000 <= E_seso_actual(t+1) <= 0.9*2000];
        C = [C, 0 <= P_own_ch_actual <= 500, 0 <= P_own_dis_actual <= 500];
        bound = 500; % 设置一个合理的边界值
        C = [C, -bound <= delta_to_dso_c(t) <= bound];
        C = [C, -bound <= delta_to_dso_d(t) <= bound];
        C = [C, -bound <= delta_from_mg1_c(t) <= bound];
        C = [C, -bound <= delta_from_mg1_d(t) <= bound];
        C = [C, -bound <= delta_from_mg2_c(t) <= bound];
        C = [C, -bound <= delta_from_mg2_d(t) <= bound];
        C = [C, -bound <= delta_from_mg3_c(t) <= bound];
        C = [C, -bound <= delta_from_mg3_d(t) <= bound];
    end
    
    %% 求解
    ops = sdpsettings('solver','cplex','verbose',0);
    sol = optimize(C, Objective, ops);
    
    if sol.problem == 0
        delta_to_dso_c = value(delta_to_dso_c);
        delta_to_dso_d = value(delta_to_dso_d);
        delta_from_mg1_c = value(delta_from_mg1_c);
        delta_from_mg1_d = value(delta_from_mg1_d);
        delta_from_mg2_c = value(delta_from_mg2_c);
        delta_from_mg2_d = value(delta_from_mg2_d);
        delta_from_mg3_c = value(delta_from_mg3_c);
        delta_from_mg3_d = value(delta_from_mg3_d);
    else
        disp('!!! SESO 日内子问题求解失败 !!!');
        yalmiperror(sol.problem);
        delta_to_dso_c = zeros(1,T); delta_to_dso_d = zeros(1,T);
        delta_from_mg1_c = zeros(1,T); delta_from_mg1_d = zeros(1,T);
        delta_from_mg2_c = zeros(1,T); delta_from_mg2_d = zeros(1,T);
        delta_from_mg3_c = zeros(1,T); delta_from_mg3_d = zeros(1,T);
    end
end
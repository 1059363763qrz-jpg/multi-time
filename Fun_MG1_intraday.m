function [delta_net_exchange, delta_lease_c, delta_lease_d] = Fun_MG1_intraday(plan_DA, mg_data_rt, k, H, ...
    delta_seso_c_k, delta_seso_d_k, lambda_c, lambda_d, rho)
%Fun_MG1_intraday: 微网1的日内滚动优化子问题 (风电)
% 版本更新: 移除了与ADMM冲突的cost_violation项

    yalmip('clear');
    T = H; 

    %% 权重系数
    w_adjust_gt   = 0.1;
    w_adjust_bat  = 0.05;
    w_adjust_dr   = 0.08;
    w_unbalance   = 100;

    %% 决策变量
    delta_e_GT    = sdpvar(1, T); delta_batc    = sdpvar(1, T);
    delta_batd    = sdpvar(1, T); delta_buy     = sdpvar(1, T);
    delta_sell    = sdpvar(1, T); delta_lease_c = sdpvar(1, T);
    delta_lease_d = sdpvar(1, T); delta_e_cut   = sdpvar(1, T);
    delta_e_tran  = sdpvar(1, T); P_unbalance   = sdpvar(1, T);

    %% 目标函数
    cost_adjust = w_adjust_gt * sum(delta_e_GT.^2) + ...
                  w_adjust_bat * sum(delta_batc.^2 + delta_batd.^2) + ...
                  w_adjust_dr * sum(delta_e_cut.^2 + delta_e_tran.^2);

    cost_unbalance = w_unbalance * sum(P_unbalance.^2);

    % =================== 关键修正：删除 cost_violation ===================
    % cost_violation 项与 ADMM 罚项冲突，必须移除
    % cost_violation = w_violation * (norm(delta_lease_c - delta_seso_c_k)^2 + ...
    %                                 norm(delta_lease_d - delta_seso_d_k)^2);
    % ==================================================================

    admm_penalty = sum(lambda_c .* (delta_lease_c - delta_seso_c_k)) ...
                 + sum(lambda_d .* (delta_lease_d - delta_seso_d_k)) ...
                 + (rho/2) * (norm(delta_lease_c - delta_seso_c_k)^2 + ...
                              norm(delta_lease_d - delta_seso_d_k)^2);
    
    Objective = cost_adjust + cost_unbalance + admm_penalty;

    %% 约束条件 (保持不变)
    C = [];
    E_bat_actual = sdpvar(1, T+1);
    hour_idx_start = floor((k-1) / 4) + 1;
    E_bat_initial_plan = 500; 
    C = [C, E_bat_actual(1) == E_bat_initial_plan]; 

    for t = 1:T
        current_k = k + t - 1;
        hour_idx = floor((current_k - 1) / 4) + 1;
        P_cut_actual = plan_DA.P_e_cut_mg1(hour_idx) + delta_e_cut(t);
        P_tran_actual = plan_DA.P_e_tran_mg1(hour_idx) + delta_e_tran(t);
        L_e_actual = mg_data_rt.L_e0_mg1(t) + P_cut_actual + P_tran_actual;
        P_gt_actual      = 0 + delta_e_GT(t);
        P_buy_actual     = max(0, -plan_DA.P_mg1_net_exchange(hour_idx)) + delta_buy(t);
        P_sell_actual    = max(0, plan_DA.P_mg1_net_exchange(hour_idx)) + delta_sell(t);
        P_batc_actual    = 0 + delta_batc(t);
        P_batd_actual    = 0 + delta_batd(t);
        P_lease_c_actual = plan_DA.P_mg1_lease_c(hour_idx) + delta_lease_c(t);
        P_lease_d_actual = plan_DA.P_mg1_lease_d(hour_idx) + delta_lease_d(t);
        total_gen = P_gt_actual + mg_data_rt.Predict_wt_mg1(t) + P_buy_actual + P_batd_actual;
        total_load = L_e_actual + P_sell_actual + P_batc_actual;
        C = [C, total_gen == total_load + P_unbalance(t)];
        eta_charge = 0.95; eta_discharge = 0.96;
        P_total_charge = P_batc_actual + P_lease_c_actual;
        P_total_discharge = P_batd_actual + P_lease_d_actual;
        C = [C, E_bat_actual(t+1) == E_bat_actual(t) + (P_total_charge * eta_charge - P_total_discharge / eta_discharge) * 0.25];
        C = [C, 100 <= E_bat_actual(t+1) <= 1500];
        C = [C, 0 <= P_total_charge <= 500, 0 <= P_total_discharge <= 500];
        C = [C, P_gt_actual >= 0, P_buy_actual >= 0, P_sell_actual >= 0];
        C = [C, P_batc_actual >= 0, P_batd_actual >= 0];
        C = [C, P_lease_c_actual >= 0, P_lease_d_actual >= 0];
        C = [C, P_cut_actual <= 0];
    end
   end_hour_idx = min(hour_idx_start + floor(H/4) - 1, 24); 
   C = [C, sum(delta_e_tran) == -sum(plan_DA.P_e_tran_mg1(hour_idx_start : end_hour_idx))];
    
    %% 求解
    ops = sdpsettings('solver','cplex','verbose',0);
    sol = optimize(C, Objective, ops);

    if sol.problem == 0
        delta_net_exchange = value(delta_sell) - value(delta_buy);
        delta_lease_c = value(delta_lease_c);
        delta_lease_d = value(delta_lease_d);
    else
        disp('!!! MG1 日内子问题求解失败 !!!');
        delta_net_exchange = zeros(1,T);
        delta_lease_c = zeros(1,T);
        delta_lease_d = zeros(1,T);
    end
end
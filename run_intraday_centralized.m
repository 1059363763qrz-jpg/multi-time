% ======================== MODIFICATION START ========================
% 增加 initial_soc_mgX 作为输入参数
function [deltas, success, diagnostics] = run_intraday_centralized(plan_DA, forecasts_RT, k, H, initial_soc_mg1, initial_soc_mg2, initial_soc_mg3, params)
% ======================== MODIFICATION END ========================
% run_intraday_centralized: 执行日内滚动调度的集中式优化版本。
%
% V3.0 更新: 增加了对日前租赁协议的违约惩罚机制。
% V4.0 更新: 增加返回t=1时刻的储能详细状态
% V5.0 更新: 修正MG储能实际充放电功率计算逻辑
% V6.0 更新: 修正MG GT实际出力计算逻辑 + 实现SOC状态继承
% V7.0 更新: 违约惩罚参数化 + 违约诊断输出

    if nargin < 8 || isempty(params)
        params = struct();
    end
    if ~isfield(params, 'violation') || isempty(params.violation)
        params.violation = 50;
    end

    yalmip('clear');
    T = H; % 优化时域长度

    fprintf('--- 执行集中式日内调度 (k=%d, H=%d, violation=%.4g) ---\n', k, T, params.violation);

    %% ======================= 1. 定义权重系数 =======================
    w = struct();
    w.dso_grid = 0.005; w.dso_dr = 0.002; w.dso_ses = 0.003;
    w.mg_gt = 0.1; w.mg_bat = 0.05; w.mg_dr = 0.08;
    w.seso_own = 0.02; w.violation = params.violation;

    %% ======================= 2. 定义所有决策变量 =======================
    % --- 调整量 ---
    delta_dso_grid = sdpvar(1, T); delta_dso_charge = sdpvar(1, T); delta_dso_discharge = sdpvar(1, T);
    delta_dso_Lshift = sdpvar(1, T); delta_dso_SIL = sdpvar(1, T);
    delta_seso_own_ch = sdpvar(1, T); delta_seso_own_dis = sdpvar(1, T);
    delta_mg1_gt = sdpvar(1, T); delta_mg1_batc = sdpvar(1, T); delta_mg1_batd = sdpvar(1, T);
    delta_mg1_buy = sdpvar(1, T); delta_mg1_sell = sdpvar(1, T); delta_mg1_lease_c = sdpvar(1, T);
    delta_mg1_lease_d = sdpvar(1, T); delta_mg1_e_cut = sdpvar(1, T); delta_mg1_e_tran = sdpvar(1, T);
    delta_mg2_gt = sdpvar(1, T); delta_mg2_batc = sdpvar(1, T); delta_mg2_batd = sdpvar(1, T);
    delta_mg2_buy = sdpvar(1, T); delta_mg2_sell = sdpvar(1, T); delta_mg2_lease_c = sdpvar(1, T);
    delta_mg2_lease_d = sdpvar(1, T); delta_mg2_e_cut = sdpvar(1, T); delta_mg2_e_tran = sdpvar(1, T);
    delta_mg3_gt = sdpvar(1, T); delta_mg3_batc = sdpvar(1, T); delta_mg3_batd = sdpvar(1, T);
    delta_mg3_buy = sdpvar(1, T); delta_mg3_sell = sdpvar(1, T); delta_mg3_lease_c = sdpvar(1, T);
    delta_mg3_lease_d = sdpvar(1, T); delta_mg3_e_cut = sdpvar(1, T); delta_mg3_e_tran = sdpvar(1, T);
    % --- 二进制 ---
    u_dso_ch = binvar(1, T); u_dso_dis = binvar(1, T); u_seso_ch = binvar(1, T); u_seso_dis = binvar(1, T);
    u_mg1_ch = binvar(1, T); u_mg1_dis = binvar(1, T); u_mg2_ch = binvar(1, T); u_mg2_dis = binvar(1, T);
    u_mg3_ch = binvar(1, T); u_mg3_dis = binvar(1, T);
    % --- 违约量 ---
    v_mg1_c = sdpvar(1, T); v_mg1_d = sdpvar(1, T); v_mg2_c = sdpvar(1, T); v_mg2_d = sdpvar(1, T);
    v_mg3_c = sdpvar(1, T); v_mg3_d = sdpvar(1, T); v_seso_c = sdpvar(1, T); v_seso_d = sdpvar(1, T);
    % --- 辅助SOC ---
    E_bat_mg1_actual = sdpvar(1, T+1); E_bat_mg2_actual = sdpvar(1, T+1);
    E_bat_mg3_actual = sdpvar(1, T+1); E_seso_actual = sdpvar(1, T+1);

    %% ======================= 3. 构建全局目标函数 =======================
    cost_adjust_dso = w.dso_grid * sum(delta_dso_grid.^2) + w.dso_dr * sum(delta_dso_Lshift.^2 + delta_dso_SIL.^2) + w.dso_ses * sum(delta_dso_charge.^2 + delta_dso_discharge.^2);
    cost_adjust_seso = w.seso_own * (sum(delta_seso_own_ch.^2) + sum(delta_seso_own_dis.^2));
    cost_adjust_mg1 = w.mg_gt * sum(delta_mg1_gt.^2) + w.mg_bat * sum(delta_mg1_batc.^2 + delta_mg1_batd.^2) + w.mg_dr * sum(delta_mg1_e_cut.^2 + delta_mg1_e_tran.^2);
    cost_adjust_mg2 = w.mg_gt * sum(delta_mg2_gt.^2) + w.mg_bat * sum(delta_mg2_batc.^2 + delta_mg2_batd.^2) + w.mg_dr * sum(delta_mg2_e_cut.^2 + delta_mg2_e_tran.^2);
    cost_adjust_mg3 = w.mg_gt * sum(delta_mg3_gt.^2) + w.mg_bat * sum(delta_mg3_batc.^2 + delta_mg3_batd.^2) + w.mg_dr * sum(delta_mg3_e_cut.^2 + delta_mg3_e_tran.^2);
    cost_adjustment = cost_adjust_dso + cost_adjust_seso + cost_adjust_mg1 + cost_adjust_mg2 + cost_adjust_mg3;
    cost_violation = w.violation * (sum(v_mg1_c) + sum(v_mg1_d) + sum(v_mg2_c) + sum(v_mg2_d) + sum(v_mg3_c) + sum(v_mg3_d) + sum(v_seso_c) + sum(v_seso_d));
    Objective = cost_adjustment + cost_violation;

    %% ======================= 4. 构建全局约束条件 =======================
    C = [];
    hour_idx_start = floor((k-1) / 4) + 1; % 日前计划中对应的起始小时索引

    % ======================== MODIFICATION START ========================
    % --- 初始状态约束 (使用传入的实际SOC) ---
    C = [C, E_bat_mg1_actual(1) == initial_soc_mg1];
    C = [C, E_bat_mg2_actual(1) == initial_soc_mg2];
    C = [C, E_bat_mg3_actual(1) == initial_soc_mg3];
    C = [C, E_seso_actual(1) == plan_DA.E_seso(hour_idx_start)]; % SESO仍然读取日前计划的SOC
    % ======================== MODIFICATION END ========================

    C = [C, v_mg1_c >= 0, v_mg1_d >= 0, v_mg2_c >= 0, v_mg2_d >= 0, v_mg3_c >= 0, v_mg3_d >= 0];
    C = [C, v_seso_c >= 0, v_seso_d >= 0];

    % --- 循环构建每个时间点的约束 ---
    P_batc_actual_mg1 = sdpvar(1, T); P_batd_actual_mg1 = sdpvar(1, T);
    P_batc_actual_mg2 = sdpvar(1, T); P_batd_actual_mg2 = sdpvar(1, T);
    P_batc_actual_mg3 = sdpvar(1, T); P_batd_actual_mg3 = sdpvar(1, T);

    for t = 1:T
        current_k = k + t - 1;
        hour_idx = floor((current_k - 1) / 4) + 1;

        % --- A. DSO 约束 ---
        P_grid_actual = plan_DA.P_dso_grid(hour_idx) + delta_dso_grid(t);
        Lshift_actual = plan_DA.Lshift_effect_dso(hour_idx) + delta_dso_Lshift(t);
        SIL_actual = plan_DA.SIL_effect_dso(hour_idx) + delta_dso_SIL(t);
        P_charge_actual_dso = plan_DA.P_dso_charge(hour_idx) + delta_dso_charge(t);
        P_discharge_actual_dso = plan_DA.P_dso_discharge(hour_idx) + delta_dso_discharge(t);
        P_mg1_net_actual = plan_DA.P_mg1_net_exchange(hour_idx) + (delta_mg1_sell(t) - delta_mg1_buy(t));
        P_mg2_net_actual = plan_DA.P_mg2_net_exchange(hour_idx) + (delta_mg2_sell(t) - delta_mg2_buy(t));
        P_mg3_net_actual = plan_DA.P_mg3_net_exchange(hour_idx) + (delta_mg3_sell(t) - delta_mg3_buy(t));
        total_gen_dso = P_grid_actual + sum(forecasts_RT.dso.P_pv(:, t)) + sum(forecasts_RT.dso.P_wt(:, t)) + P_discharge_actual_dso - P_mg1_net_actual - P_mg2_net_actual - P_mg3_net_actual;
        total_load_dso = sum(forecasts_RT.dso.PLOAD(:, t)) + Lshift_actual - SIL_actual + P_charge_actual_dso;
        C = [C, total_gen_dso == total_load_dso + plan_DA.P_loss(hour_idx)];
        P_dso_seso_max = 1000;
        C = [C, 0 <= P_charge_actual_dso <= P_dso_seso_max * u_dso_ch(t)];
        C = [C, 0 <= P_discharge_actual_dso <= P_dso_seso_max * u_dso_dis(t)];
        C = [C, u_dso_ch(t) + u_dso_dis(t) <= 1];

        % --- B. SESO 约束 ---
        P_seso_ch_actual = plan_DA.P_seso_ch(hour_idx) + delta_seso_own_ch(t);
        P_seso_dis_actual = plan_DA.P_seso_dis(hour_idx) + delta_seso_own_dis(t);
        C = [C, E_seso_actual(t+1) == E_seso_actual(t) + (P_seso_ch_actual*0.95 - P_seso_dis_actual/0.95)*0.25];
        C = [C, 0.1*2000 <= E_seso_actual(t+1) <= 0.9*2000];
        P_seso_max = 5000;
        C = [C, 0 <= P_seso_ch_actual <= P_seso_max * u_seso_ch(t)];
        C = [C, 0 <= P_seso_dis_actual <= P_seso_max * u_seso_dis(t)];
        C = [C, u_seso_ch(t) + u_seso_dis(t) <= 1];

        % --- SESO 功率池与违约约束 ---
        total_charge_collected_actual = (plan_DA.P_mg1_lease_c(hour_idx) + delta_mg1_lease_c(t)) + (plan_DA.P_mg2_lease_c(hour_idx) + delta_mg2_lease_c(t)) + (plan_DA.P_mg3_lease_c(hour_idx) + delta_mg3_lease_c(t)) + P_seso_ch_actual;
        total_discharge_collected_actual = (plan_DA.P_mg1_lease_d(hour_idx) + delta_mg1_lease_d(t)) + (plan_DA.P_mg2_lease_d(hour_idx) + delta_mg2_lease_d(t)) + (plan_DA.P_mg3_lease_d(hour_idx) + delta_mg3_lease_d(t)) + P_seso_dis_actual;
        C = [C, P_charge_actual_dso == total_charge_collected_actual];
        C = [C, P_discharge_actual_dso == total_discharge_collected_actual];
        C = [C, (plan_DA.P_mg1_lease_c(hour_idx) + delta_mg1_lease_c(t)) + v_mg1_c(t) >= plan_DA.P_mg1_lease_c(hour_idx)];
        C = [C, (plan_DA.P_mg1_lease_d(hour_idx) + delta_mg1_lease_d(t)) + v_mg1_d(t) >= plan_DA.P_mg1_lease_d(hour_idx)];
        C = [C, (plan_DA.P_mg2_lease_c(hour_idx) + delta_mg2_lease_c(t)) + v_mg2_c(t) >= plan_DA.P_mg2_lease_c(hour_idx)];
        C = [C, (plan_DA.P_mg2_lease_d(hour_idx) + delta_mg2_lease_d(t)) + v_mg2_d(t) >= plan_DA.P_mg2_lease_d(hour_idx)];
        C = [C, (plan_DA.P_mg3_lease_c(hour_idx) + delta_mg3_lease_c(t)) + v_mg3_c(t) >= plan_DA.P_mg3_lease_c(hour_idx)];
        C = [C, (plan_DA.P_mg3_lease_d(hour_idx) + delta_mg3_lease_d(t)) + v_mg3_d(t) >= plan_DA.P_mg3_lease_d(hour_idx)];
        C = [C, P_seso_ch_actual + v_seso_c(t) >= plan_DA.P_seso_ch(hour_idx)];
        C = [C, P_seso_dis_actual + v_seso_d(t) >= plan_DA.P_seso_dis(hour_idx)];

        % --- C, D, E. 微网约束 ---
        [C, P_batc_actual_mg1(t), P_batd_actual_mg1(t)] = model_mg_constraints(C, t, k, plan_DA, forecasts_RT.mg, 1, delta_mg1_gt, delta_mg1_batc, delta_mg1_batd, delta_mg1_buy, delta_mg1_sell, delta_mg1_lease_c, delta_mg1_lease_d, delta_mg1_e_cut, delta_mg1_e_tran, E_bat_mg1_actual, u_mg1_ch, u_mg1_dis);
        [C, P_batc_actual_mg2(t), P_batd_actual_mg2(t)] = model_mg_constraints(C, t, k, plan_DA, forecasts_RT.mg, 2, delta_mg2_gt, delta_mg2_batc, delta_mg2_batd, delta_mg2_buy, delta_mg2_sell, delta_mg2_lease_c, delta_mg2_lease_d, delta_mg2_e_cut, delta_mg2_e_tran, E_bat_mg2_actual, u_mg2_ch, u_mg2_dis);
        [C, P_batc_actual_mg3(t), P_batd_actual_mg3(t)] = model_mg_constraints(C, t, k, plan_DA, forecasts_RT.mg, 3, delta_mg3_gt, delta_mg3_batc, delta_mg3_batd, delta_mg3_buy, delta_mg3_sell, delta_mg3_lease_c, delta_mg3_lease_d, delta_mg3_e_cut, delta_mg3_e_tran, E_bat_mg3_actual, u_mg3_ch, u_mg3_dis);
    end

    %% ======================= 5. 求解与输出 =======================
    diagnostics = struct();
    diagnostics.params = params;
    diagnostics.success = false;
    diagnostics.violation_series = struct('v_mg1_c', [], 'v_mg1_d', [], 'v_mg2_c', [], 'v_mg2_d', [], 'v_mg3_c', [], 'v_mg3_d', [], 'v_seso_c', [], 'v_seso_d', []);
    diagnostics.cost = struct('adjustment_cost', NaN, 'penalty_cost', NaN, 'objective_total', NaN, 'actual_operating_cost_excl_penalty', NaN);

    ops = sdpsettings('solver', 'gurobi', 'verbose', 0);
    tic;
    sol = optimize(C, Objective, ops);
    diagnostics.solve_time = toc;

    if sol.problem == 0
        fprintf('    集中式调度求解成功！\n');
        success = true;
        diagnostics.success = true;

        deltas.dso_grid = value(delta_dso_grid);
        deltas.dso_seso_c = value(delta_dso_charge);
        deltas.dso_seso_d = value(delta_dso_discharge);

        deltas.mg1_net = value(delta_mg1_sell) - value(delta_mg1_buy);
        deltas.mg1_seso_c = value(delta_mg1_lease_c);
        deltas.mg1_seso_d = value(delta_mg1_lease_d);
        deltas.mg1_gt_delta = value(delta_mg1_gt);
        deltas.mg1_e_cut_delta = value(delta_mg1_e_cut);
        deltas.mg1_e_tran_delta = value(delta_mg1_e_tran);

        deltas.mg2_net = value(delta_mg2_sell) - value(delta_mg2_buy);
        deltas.mg2_seso_c = value(delta_mg2_lease_c);
        deltas.mg2_seso_d = value(delta_mg2_lease_d);
        deltas.mg2_gt_delta = value(delta_mg2_gt);
        deltas.mg2_e_cut_delta = value(delta_mg2_e_cut);
        deltas.mg2_e_tran_delta = value(delta_mg2_e_tran);

        deltas.mg3_net = value(delta_mg3_sell) - value(delta_mg3_buy);
        deltas.mg3_seso_c = value(delta_mg3_lease_c);
        deltas.mg3_seso_d = value(delta_mg3_lease_d);
        deltas.mg3_gt_delta = value(delta_mg3_gt);
        deltas.mg3_e_cut_delta = value(delta_mg3_e_cut);
        deltas.mg3_e_tran_delta = value(delta_mg3_e_tran);

        % 返回t=1时刻的储能详细状态
        hour_idx_t1 = floor(((k + 1 - 1) - 1) / 4) + 1;
        deltas.E_seso_actual_t1 = value(E_seso_actual(2));
        deltas.P_seso_ch_actual_t1 = value(plan_DA.P_seso_ch(hour_idx_t1) + delta_seso_own_ch(1));
        deltas.P_seso_dis_actual_t1 = value(plan_DA.P_seso_dis(hour_idx_t1) + delta_seso_own_dis(1));
        deltas.E_mg1_actual_t1 = value(E_bat_mg1_actual(2));
        deltas.P_mg1_batc_actual_t1 = value(P_batc_actual_mg1(1));
        deltas.P_mg1_batd_actual_t1 = value(P_batd_actual_mg1(1));
        deltas.E_mg2_actual_t1 = value(E_bat_mg2_actual(2));
        deltas.P_mg2_batc_actual_t1 = value(P_batc_actual_mg2(1));
        deltas.P_mg2_batd_actual_t1 = value(P_batd_actual_mg2(1));
        deltas.E_mg3_actual_t1 = value(E_bat_mg3_actual(2));
        deltas.P_mg3_batc_actual_t1 = value(P_batc_actual_mg3(1));
        deltas.P_mg3_batd_actual_t1 = value(P_batd_actual_mg3(1));

        % 返回t=1时刻关键可调量
        h1 = floor((k - 1) / 4) + 1;
        deltas.P_mg1_gt_actual_t1 = plan_DA.P_mg1_gt(h1) + deltas.mg1_gt_delta(1);
        deltas.P_mg2_gt_actual_t1 = plan_DA.P_mg2_gt(h1) + deltas.mg2_gt_delta(1);
        deltas.P_mg3_gt_actual_t1 = plan_DA.P_mg3_gt(h1) + deltas.mg3_gt_delta(1);
        deltas.P_mg1_e_cut_actual_t1 = plan_DA.P_e_cut_mg1(h1) + deltas.mg1_e_cut_delta(1);
        deltas.P_mg2_e_cut_actual_t1 = plan_DA.P_e_cut_mg2(h1) + deltas.mg2_e_cut_delta(1);
        deltas.P_mg3_e_cut_actual_t1 = plan_DA.P_e_cut_mg3(h1) + deltas.mg3_e_cut_delta(1);
        deltas.P_mg1_e_tran_actual_t1 = plan_DA.P_e_tran_mg1(h1) + deltas.mg1_e_tran_delta(1);
        deltas.P_mg2_e_tran_actual_t1 = plan_DA.P_e_tran_mg2(h1) + deltas.mg2_e_tran_delta(1);
        deltas.P_mg3_e_tran_actual_t1 = plan_DA.P_e_tran_mg3(h1) + deltas.mg3_e_tran_delta(1);

        diagnostics.violation_series.v_mg1_c = value(v_mg1_c);
        diagnostics.violation_series.v_mg1_d = value(v_mg1_d);
        diagnostics.violation_series.v_mg2_c = value(v_mg2_c);
        diagnostics.violation_series.v_mg2_d = value(v_mg2_d);
        diagnostics.violation_series.v_mg3_c = value(v_mg3_c);
        diagnostics.violation_series.v_mg3_d = value(v_mg3_d);
        diagnostics.violation_series.v_seso_c = value(v_seso_c);
        diagnostics.violation_series.v_seso_d = value(v_seso_d);

        diagnostics.cost.adjustment_cost = value(cost_adjustment);
        diagnostics.cost.penalty_cost = value(cost_violation);
        diagnostics.cost.objective_total = value(Objective);
        diagnostics.cost.actual_operating_cost_excl_penalty = diagnostics.cost.adjustment_cost;
    else
        fprintf('!!! 警告: 集中式日内调度求解失败 !!!\n');
        disp(yalmiperror(sol.problem));
        success = false;
        deltas = [];
        diagnostics.error_code = sol.problem;
        diagnostics.error_message = yalmiperror(sol.problem);
    end
end

%% 辅助函数 (V4: 修正 P_gt_actual 计算)
function [C, P_batc_actual, P_batd_actual] = model_mg_constraints(C, t, k, plan_DA, mg_rt, mg_idx, ...
    delta_gt, delta_batc, delta_batd, delta_buy, delta_sell, ...
    delta_lease_c, delta_lease_d, delta_e_cut, delta_e_tran, E_bat_actual, ...
    u_mg_ch, u_mg_dis)

    hour_idx = floor(((k+t-1)-1) / 4) + 1;

    P_e_cut_plan = plan_DA.(['P_e_cut_mg', num2str(mg_idx)])(hour_idx);
    P_e_tran_plan = plan_DA.(['P_e_tran_mg', num2str(mg_idx)])(hour_idx);
    L_e0_rt = mg_rt.(['L_e0_mg', num2str(mg_idx)])(t);
    P_net_ex_plan = plan_DA.(['P_mg', num2str(mg_idx), '_net_exchange'])(hour_idx);
    P_lease_c_plan = plan_DA.(['P_mg', num2str(mg_idx), '_lease_c'])(hour_idx);
    P_lease_d_plan = plan_DA.(['P_mg', num2str(mg_idx), '_lease_d'])(hour_idx);
    P_batc_plan = plan_DA.(['P_mg', num2str(mg_idx), '_batc'])(hour_idx);
    P_batd_plan = plan_DA.(['P_mg', num2str(mg_idx), '_batd'])(hour_idx);

    % ======================== MODIFICATION START ========================
    % 从plan_DA中读取日前计划的 P_gt
    P_gt_plan = plan_DA.(['P_mg', num2str(mg_idx), '_gt'])(hour_idx);
    % ======================== MODIFICATION END ========================

    if mg_idx == 1, P_re_rt = mg_rt.Predict_wt_mg1(t);
    else, P_re_rt = mg_rt.(['Predict_pv_mg', num2str(mg_idx)])(t); end

    P_cut_actual = P_e_cut_plan + delta_e_cut(t);
    P_tran_actual = P_e_tran_plan + delta_e_tran(t);
    L_e_actual = L_e0_rt + P_cut_actual + P_tran_actual;
    % ======================== MODIFICATION START ========================
    P_gt_actual = P_gt_plan + delta_gt(t); % <--- 使用日前计划值 + 调整量
    % ======================== MODIFICATION END ========================
    P_buy_actual = max(0, -P_net_ex_plan) + delta_buy(t);
    P_sell_actual = max(0, P_net_ex_plan) + delta_sell(t);
    P_batc_actual = P_batc_plan + delta_batc(t);
    P_batd_actual = P_batd_plan + delta_batd(t);
    P_lease_c_actual = P_lease_c_plan + delta_lease_c(t);
    P_lease_d_actual = P_lease_d_plan + delta_lease_d(t);

    total_gen = P_gt_actual + P_re_rt + P_buy_actual + P_batd_actual;
    total_load = L_e_actual + P_sell_actual + P_batc_actual;
    C = [C, total_gen == total_load];

    P_total_charge = P_batc_actual + P_lease_c_actual;
    P_total_discharge = P_batd_actual + P_lease_d_actual;
    C = [C, E_bat_actual(t+1) == E_bat_actual(t) + (P_total_charge * 0.95 - P_total_discharge / 0.96) * 0.25];
    C = [C, 500 <= E_bat_actual(t+1) <= 1800]; % MG储能容量约束

    P_mg_bat_max = 500; % MG储能最大功率
    C = [C, 0 <= P_total_charge <= P_mg_bat_max * u_mg_ch(t)];
    C = [C, 0 <= P_total_discharge <= P_mg_bat_max * u_mg_dis(t)];
    C = [C, u_mg_ch(t) + u_mg_dis(t) <= 1]; % 充放电互斥

    C = [C, P_gt_actual >= 0, P_buy_actual >= 0, P_sell_actual >= 0];
    C = [C, P_batc_actual >= 0, P_batd_actual >= 0];
    C = [C, P_lease_c_actual >= 0, P_lease_d_actual >= 0];
    C = [C, P_cut_actual <= 0];
end

function [actual_operation, run_stats] = run_intraday_process(plan_DA, dso_data, mg_data, T_RT, H, params)
% run_intraday_process: 执行日内滚动调度公共流程。
% 说明：仅拆分流程，不修改优化模型与约束。

    if nargin < 6 || isempty(params)
        params = struct();
    end
    if ~isfield(params, 'violation') || isempty(params.violation)
        params.violation = 50;
    end
    if ~isfield(params, 'fail_fast') || isempty(params.fail_fast)
        params.fail_fast = false;
    end

    n_rt = length(T_RT);

    % 初始化用于存储日内实际运行结果的变量
    actual_operation = struct();
    actual_operation.P_dso_grid = zeros(1, n_rt);
    actual_operation.P_dso_seso_c = zeros(1, n_rt);
    actual_operation.P_dso_seso_d = zeros(1, n_rt);
    actual_operation.P_mg1_net = zeros(1, n_rt);
    actual_operation.P_mg2_net = zeros(1, n_rt);
    actual_operation.P_mg3_net = zeros(1, n_rt);
    actual_operation.P_mg1_seso_c = zeros(1, n_rt);
    actual_operation.P_mg1_seso_d = zeros(1, n_rt);
    actual_operation.P_mg2_seso_c = zeros(1, n_rt);
    actual_operation.P_mg2_seso_d = zeros(1, n_rt);
    actual_operation.P_mg3_seso_c = zeros(1, n_rt);
    actual_operation.P_mg3_seso_d = zeros(1, n_rt);
    actual_operation.E_seso = zeros(1, n_rt);
    actual_operation.P_seso_ch = zeros(1, n_rt);
    actual_operation.P_seso_dis = zeros(1, n_rt);
    actual_operation.E_mg1 = zeros(1, n_rt);
    actual_operation.P_mg1_batc = zeros(1, n_rt);
    actual_operation.P_mg1_batd = zeros(1, n_rt);
    actual_operation.E_mg2 = zeros(1, n_rt);
    actual_operation.P_mg2_batc = zeros(1, n_rt);
    actual_operation.P_mg2_batd = zeros(1, n_rt);
    actual_operation.E_mg3 = zeros(1, n_rt);
    actual_operation.P_mg3_batc = zeros(1, n_rt);
    actual_operation.P_mg3_batd = zeros(1, n_rt);

    % 增加分析时序
    actual_operation.P_mg1_gt = zeros(1, n_rt);
    actual_operation.P_mg2_gt = zeros(1, n_rt);
    actual_operation.P_mg3_gt = zeros(1, n_rt);
    actual_operation.P_mg1_e_cut = zeros(1, n_rt);
    actual_operation.P_mg2_e_cut = zeros(1, n_rt);
    actual_operation.P_mg3_e_cut = zeros(1, n_rt);
    actual_operation.P_mg1_e_tran = zeros(1, n_rt);
    actual_operation.P_mg2_e_tran = zeros(1, n_rt);
    actual_operation.P_mg3_e_tran = zeros(1, n_rt);

    run_stats = struct();
    run_stats.params = params;
    run_stats.step_success = false(1, n_rt);
    run_stats.failed_steps_count = 0;
    run_stats.solve_time_total = 0;
    run_stats.solve_time_series = zeros(1, n_rt);
    run_stats.cost_step = struct('adjustment_cost', zeros(1, n_rt), 'penalty_cost', zeros(1, n_rt), 'objective_total', zeros(1, n_rt), 'actual_operating_cost_excl_penalty', zeros(1, n_rt));
    run_stats.violation_series = struct('v_mg1_c', zeros(1, n_rt), 'v_mg1_d', zeros(1, n_rt), ...
                                        'v_mg2_c', zeros(1, n_rt), 'v_mg2_d', zeros(1, n_rt), ...
                                        'v_mg3_c', zeros(1, n_rt), 'v_mg3_d', zeros(1, n_rt), ...
                                        'v_seso_c', zeros(1, n_rt), 'v_seso_d', zeros(1, n_rt));

    % 为 k=1 设置初始SOC
    initial_soc_mg1 = 800;
    initial_soc_mg2 = 800;
    initial_soc_mg3 = 800;

    % 模拟日内滚动过程
    for k = 1:n_rt
        fprintf('--- 正在处理日内时间点 k = %d ---\n', k);

        % a. 确定当前滚动窗口
        window_idx = k:min(k + H - 1, n_rt);
        window_len = length(window_idx);

        % b. 获取当前窗口的实时预测数据
        current_forecasts_RT.dso = filter_data_by_window(dso_data.forecasts_RT, window_idx, n_rt);
        current_forecasts_RT.mg = filter_data_by_window(mg_data.forecasts_RT, window_idx, n_rt);

        % c. 获取上一个时间点的实际SOC作为当前优化的初始SOC
        if k == 1
            current_initial_soc_mg1 = initial_soc_mg1;
            current_initial_soc_mg2 = initial_soc_mg2;
            current_initial_soc_mg3 = initial_soc_mg3;
        else
            current_initial_soc_mg1 = actual_operation.E_mg1(k-1);
            current_initial_soc_mg2 = actual_operation.E_mg2(k-1);
            current_initial_soc_mg3 = actual_operation.E_mg3(k-1);
        end

        % d. 执行日内滚动优化 (传入初始SOC)
        [deltas, success, diagnostics] = run_intraday_centralized(plan_DA, current_forecasts_RT, k, window_len, ...
                                                                  current_initial_soc_mg1, current_initial_soc_mg2, current_initial_soc_mg3, params);

        run_stats.step_success(k) = success;
        if isfield(diagnostics, 'solve_time')
            run_stats.solve_time_series(k) = diagnostics.solve_time;
            run_stats.solve_time_total = run_stats.solve_time_total + diagnostics.solve_time;
        end

        % 获取当前时间点对应的日前小时索引
        hour_idx_DA = floor((k-1) / 4) + 1;

        if ~success || isempty(deltas)
            run_stats.failed_steps_count = run_stats.failed_steps_count + 1;
            msg = sprintf('在时间点 k = %d 的日内调度求解失败。', k);
            if params.fail_fast
                error('%s 为避免NaN传播，已fail-fast终止。', msg);
            else
                warning('%s 将使用日前计划值 (调整量为0)。', msg);
            end
            deltas = get_zero_deltas(window_len);

            % 强制使用日前计划值填充储能状态
            actual_operation.E_seso(k) = plan_DA.E_seso(hour_idx_DA);
            actual_operation.P_seso_ch(k) = plan_DA.P_seso_ch(hour_idx_DA);
            actual_operation.P_seso_dis(k) = plan_DA.P_seso_dis(hour_idx_DA);

            % 对于MG储能，如果优化失败，SOC保持上一个时刻的值
            if k == 1
                actual_operation.E_mg1(k) = initial_soc_mg1;
                actual_operation.E_mg2(k) = initial_soc_mg2;
                actual_operation.E_mg3(k) = initial_soc_mg3;
            else
                actual_operation.E_mg1(k) = actual_operation.E_mg1(k-1);
                actual_operation.E_mg2(k) = actual_operation.E_mg2(k-1);
                actual_operation.E_mg3(k) = actual_operation.E_mg3(k-1);
            end
            actual_operation.P_mg1_batc(k) = plan_DA.P_mg1_batc(hour_idx_DA);
            actual_operation.P_mg1_batd(k) = plan_DA.P_mg1_batd(hour_idx_DA);
            actual_operation.P_mg2_batc(k) = plan_DA.P_mg2_batc(hour_idx_DA);
            actual_operation.P_mg2_batd(k) = plan_DA.P_mg2_batd(hour_idx_DA);
            actual_operation.P_mg3_batc(k) = plan_DA.P_mg3_batc(hour_idx_DA);
            actual_operation.P_mg3_batd(k) = plan_DA.P_mg3_batd(hour_idx_DA);

            actual_operation.P_mg1_gt(k) = plan_DA.P_mg1_gt(hour_idx_DA);
            actual_operation.P_mg2_gt(k) = plan_DA.P_mg2_gt(hour_idx_DA);
            actual_operation.P_mg3_gt(k) = plan_DA.P_mg3_gt(hour_idx_DA);
            actual_operation.P_mg1_e_cut(k) = plan_DA.P_e_cut_mg1(hour_idx_DA);
            actual_operation.P_mg2_e_cut(k) = plan_DA.P_e_cut_mg2(hour_idx_DA);
            actual_operation.P_mg3_e_cut(k) = plan_DA.P_e_cut_mg3(hour_idx_DA);
            actual_operation.P_mg1_e_tran(k) = plan_DA.P_e_tran_mg1(hour_idx_DA);
            actual_operation.P_mg2_e_tran(k) = plan_DA.P_e_tran_mg2(hour_idx_DA);
            actual_operation.P_mg3_e_tran(k) = plan_DA.P_e_tran_mg3(hour_idx_DA);
        else
            run_stats.cost_step.adjustment_cost(k) = diagnostics.cost.adjustment_cost;
            run_stats.cost_step.penalty_cost(k) = diagnostics.cost.penalty_cost;
            run_stats.cost_step.objective_total(k) = diagnostics.cost.objective_total;
            run_stats.cost_step.actual_operating_cost_excl_penalty(k) = diagnostics.cost.actual_operating_cost_excl_penalty;
            run_stats.violation_series.v_mg1_c(k) = diagnostics.violation_series.v_mg1_c(1);
            run_stats.violation_series.v_mg1_d(k) = diagnostics.violation_series.v_mg1_d(1);
            run_stats.violation_series.v_mg2_c(k) = diagnostics.violation_series.v_mg2_c(1);
            run_stats.violation_series.v_mg2_d(k) = diagnostics.violation_series.v_mg2_d(1);
            run_stats.violation_series.v_mg3_c(k) = diagnostics.violation_series.v_mg3_c(1);
            run_stats.violation_series.v_mg3_d(k) = diagnostics.violation_series.v_mg3_d(1);
            run_stats.violation_series.v_seso_c(k) = diagnostics.violation_series.v_seso_c(1);
            run_stats.violation_series.v_seso_d(k) = diagnostics.violation_series.v_seso_d(1);
        end

        % e. 记录第一个时间点的实际执行结果
        actual_operation.P_dso_grid(k) = plan_DA.P_dso_grid(hour_idx_DA) + deltas.dso_grid(1);
        actual_operation.P_dso_seso_c(k) = plan_DA.P_dso_charge(hour_idx_DA) + deltas.dso_seso_c(1);
        actual_operation.P_dso_seso_d(k) = plan_DA.P_dso_discharge(hour_idx_DA) + deltas.dso_seso_d(1);

        actual_operation.P_mg1_net(k) = plan_DA.P_mg1_net_exchange(hour_idx_DA) + deltas.mg1_net(1);
        actual_operation.P_mg1_seso_c(k) = plan_DA.P_mg1_lease_c(hour_idx_DA) + deltas.mg1_seso_c(1);
        actual_operation.P_mg1_seso_d(k) = plan_DA.P_mg1_lease_d(hour_idx_DA) + deltas.mg1_seso_d(1);

        actual_operation.P_mg2_net(k) = plan_DA.P_mg2_net_exchange(hour_idx_DA) + deltas.mg2_net(1);
        actual_operation.P_mg2_seso_c(k) = plan_DA.P_mg2_lease_c(hour_idx_DA) + deltas.mg2_seso_c(1);
        actual_operation.P_mg2_seso_d(k) = plan_DA.P_mg2_lease_d(hour_idx_DA) + deltas.mg2_seso_d(1);

        actual_operation.P_mg3_net(k) = plan_DA.P_mg3_net_exchange(hour_idx_DA) + deltas.mg3_net(1);
        actual_operation.P_mg3_seso_c(k) = plan_DA.P_mg3_lease_c(hour_idx_DA) + deltas.mg3_seso_c(1);
        actual_operation.P_mg3_seso_d(k) = plan_DA.P_mg3_lease_d(hour_idx_DA) + deltas.mg3_seso_d(1);

        % 记录详细的储能状态 (仅在成功时执行)
        if success
            actual_operation.E_seso(k) = deltas.E_seso_actual_t1;
            actual_operation.P_seso_ch(k) = deltas.P_seso_ch_actual_t1;
            actual_operation.P_seso_dis(k) = deltas.P_seso_dis_actual_t1;

            actual_operation.E_mg1(k) = deltas.E_mg1_actual_t1;
            actual_operation.P_mg1_batc(k) = deltas.P_mg1_batc_actual_t1;
            actual_operation.P_mg1_batd(k) = deltas.P_mg1_batd_actual_t1;

            actual_operation.E_mg2(k) = deltas.E_mg2_actual_t1;
            actual_operation.P_mg2_batc(k) = deltas.P_mg2_batc_actual_t1;
            actual_operation.P_mg2_batd(k) = deltas.P_mg2_batd_actual_t1;

            actual_operation.E_mg3(k) = deltas.E_mg3_actual_t1;
            actual_operation.P_mg3_batc(k) = deltas.P_mg3_batc_actual_t1;
            actual_operation.P_mg3_batd(k) = deltas.P_mg3_batd_actual_t1;

            if isfield(deltas, 'P_mg1_gt_actual_t1')
                actual_operation.P_mg1_gt(k) = deltas.P_mg1_gt_actual_t1;
                actual_operation.P_mg2_gt(k) = deltas.P_mg2_gt_actual_t1;
                actual_operation.P_mg3_gt(k) = deltas.P_mg3_gt_actual_t1;
                actual_operation.P_mg1_e_cut(k) = deltas.P_mg1_e_cut_actual_t1;
                actual_operation.P_mg2_e_cut(k) = deltas.P_mg2_e_cut_actual_t1;
                actual_operation.P_mg3_e_cut(k) = deltas.P_mg3_e_cut_actual_t1;
                actual_operation.P_mg1_e_tran(k) = deltas.P_mg1_e_tran_actual_t1;
                actual_operation.P_mg2_e_tran(k) = deltas.P_mg2_e_tran_actual_t1;
                actual_operation.P_mg3_e_tran(k) = deltas.P_mg3_e_tran_actual_t1;
            end
        end

        fprintf('时间点 k = %d 处理完毕。\n', k);
    end

    run_stats.success = all(run_stats.step_success);
    if n_rt > 0
        run_stats.solve_time_avg = run_stats.solve_time_total / n_rt;
    else
        run_stats.solve_time_avg = 0;
    end
end

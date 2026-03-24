function actual_operation = run_intraday_process(plan_DA, dso_data, mg_data, T_RT, H)
% run_intraday_process: 执行日内滚动调度公共流程。
% 说明：仅拆分流程，不修改优化模型与约束。

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
        [deltas, success] = run_intraday_centralized(plan_DA, current_forecasts_RT, k, window_len, ...
                                                     current_initial_soc_mg1, current_initial_soc_mg2, current_initial_soc_mg3);

        % 获取当前时间点对应的日前小时索引
        hour_idx_DA = floor((k-1) / 4) + 1;

        if ~success || isempty(deltas)
            warning('在时间点 k = %d 的日内调度求解失败！将使用日前计划值 (调整量为0)。', k);
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
        end

        fprintf('时间点 k = %d 处理完毕。\n', k);
    end
end

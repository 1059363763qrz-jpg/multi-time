%% 主函数：计及共享储能的配电网-多微网多时间尺度协同调度 (修正版)
clc; clear; close all;

fprintf('======================================================\n');
fprintf('      配电网-微网多时间尺度协同调度程序 V2.0 (集中式日内) \n');
fprintf('======================================================\n\n');

%% 1. 系统参数定义
T_DA = 1:24;
T_RT = (15/60):(15/60):24; % 修正时间轴以匹配96个点
H = 4;                     % 日内滚动优化窗口长度 (4个15分钟，即1小时)

%% 2. 数据准备
fprintf('(1/4) 正在准备日前及日内数据...\n');
dso_data = prepare_data_DSO(T_DA, T_RT);
mg_data = prepare_data_MGs(T_DA, T_RT);
fprintf('数据准备完毕。\n\n');

%% 3. 日前调度 (Day-Ahead Scheduling)
fprintf('(2/4) 正在执行日前经济调度...\n');
[plan_DA, admm_results_DA] = run_day_ahead_scheduling(dso_data.forecasts_DA, mg_data.forecasts_DA);
if isempty(plan_DA)
    error('日前调度失败，程序终止。');
end
fprintf('日前调度完成，已生成24小时运行计划。\n\n');

%% 4. 日内滚动调度 (Intra-day Rolling Scheduling)
fprintf('(3/4) 开始执行日内滚动调度 (总共 %d 个时间点)...\n', length(T_RT));

% 初始化用于存储日内实际运行结果的变量
actual_operation = struct();
actual_operation.P_dso_grid = zeros(1, 96);
actual_operation.P_dso_seso_c = zeros(1, 96);
actual_operation.P_dso_seso_d = zeros(1, 96);
actual_operation.P_mg1_net = zeros(1, 96);
actual_operation.P_mg2_net = zeros(1, 96);
actual_operation.P_mg3_net = zeros(1, 96);
actual_operation.P_mg1_seso_c = zeros(1, 96);
actual_operation.P_mg1_seso_d = zeros(1, 96);
actual_operation.P_mg2_seso_c = zeros(1, 96);
actual_operation.P_mg2_seso_d = zeros(1, 96);
actual_operation.P_mg3_seso_c = zeros(1, 96);
actual_operation.P_mg3_seso_d = zeros(1, 96);
actual_operation.E_seso = zeros(1, 96);
actual_operation.P_seso_ch = zeros(1, 96);
actual_operation.P_seso_dis = zeros(1, 96);
actual_operation.E_mg1 = zeros(1, 96);
actual_operation.P_mg1_batc = zeros(1, 96);
actual_operation.P_mg1_batd = zeros(1, 96);
actual_operation.E_mg2 = zeros(1, 96);
actual_operation.P_mg2_batc = zeros(1, 96);
actual_operation.P_mg2_batd = zeros(1, 96);
actual_operation.E_mg3 = zeros(1, 96);
actual_operation.P_mg3_batc = zeros(1, 96);
actual_operation.P_mg3_batd = zeros(1, 96);

% ======================== MODIFICATION START ========================
% 为 k=1 设置初始SOC (假设与日前计划一致或固定值)
initial_soc_mg1 = 800; % Or plan_DA.E_mg1(1);
initial_soc_mg2 = 800; % Or plan_DA.E_mg2(1);
initial_soc_mg3 = 800; % Or plan_DA.E_mg3(1);
% 注意：SESO的初始SOC由run_intraday_centralized内部从plan_DA读取，无需在此设置
% ======================== MODIFICATION END ========================


% 模拟日内滚动过程
for k = 1:96
    fprintf('--- 正在处理日内时间点 k = %d ---\n', k);
    
    % a. 确定当前滚动窗口
    window_idx = k:min(k + H - 1, 96);
    window_len = length(window_idx);

    % b. 获取当前窗口的实时预测数据
    current_forecasts_RT.dso = filter_data_by_window(dso_data.forecasts_RT, window_idx);
    current_forecasts_RT.mg  = filter_data_by_window(mg_data.forecasts_RT, window_idx);

    % ======================== MODIFICATION START ========================
    % c. 获取上一个时间点的实际SOC作为当前优化的初始SOC
    if k == 1
        current_initial_soc_mg1 = initial_soc_mg1;
        current_initial_soc_mg2 = initial_soc_mg2;
        current_initial_soc_mg3 = initial_soc_mg3;
    else
        % 从上一个时间点k-1的实际运行结果中获取SOC
        current_initial_soc_mg1 = actual_operation.E_mg1(k-1);
        current_initial_soc_mg2 = actual_operation.E_mg2(k-1);
        current_initial_soc_mg3 = actual_operation.E_mg3(k-1);
    end
    
    % d. 执行日内滚动优化 (传入初始SOC)
    [deltas, success] = run_intraday_centralized(plan_DA, current_forecasts_RT, k, window_len, ...
                                                 current_initial_soc_mg1, current_initial_soc_mg2, current_initial_soc_mg3);
    % ======================== MODIFICATION END ========================
                                                 
    % 获取当前时间点对应的日前小时索引
    hour_idx_DA = floor((k-1) / 4) + 1;
    
    if ~success || isempty(deltas)
        warning('在时间点 k = %d 的日内调度求解失败！将使用日前计划值 (调整量为0)。', k);
        deltas = get_zero_deltas(window_len);
        
        % 强制使用日前计划值填充储能状态 (注意：SOC需要特殊处理，因为是状态量)
        actual_operation.E_seso(k) = plan_DA.E_seso(hour_idx_DA); % 简化处理，实际应模拟演变
        actual_operation.P_seso_ch(k) = plan_DA.P_seso_ch(hour_idx_DA);
        actual_operation.P_seso_dis(k) = plan_DA.P_seso_dis(hour_idx_DA);
        
        % 对于MG储能，如果优化失败，SOC保持上一个时刻的值
        if k==1
             actual_operation.E_mg1(k) = initial_soc_mg1;
             actual_operation.E_mg2(k) = initial_soc_mg2;
             actual_operation.E_mg3(k) = initial_soc_mg3;
        else
             actual_operation.E_mg1(k) = actual_operation.E_mg1(k-1);
             actual_operation.E_mg2(k) = actual_operation.E_mg2(k-1);
             actual_operation.E_mg3(k) = actual_operation.E_mg3(k-1);
        end
        actual_operation.P_mg1_batc(k) = plan_DA.P_mg1_batc(hour_idx_DA); % 功率使用日前计划
        actual_operation.P_mg1_batd(k) = plan_DA.P_mg1_batd(hour_idx_DA);
        actual_operation.P_mg2_batc(k) = plan_DA.P_mg2_batc(hour_idx_DA);
        actual_operation.P_mg2_batd(k) = plan_DA.P_mg2_batd(hour_idx_DA);
        actual_operation.P_mg3_batc(k) = plan_DA.P_mg3_batc(hour_idx_DA);
        actual_operation.P_mg3_batd(k) = plan_DA.P_mg3_batd(hour_idx_DA);

    end

    % e. 记录第一个时间点的实际执行结果
    actual_operation.P_dso_grid(k)     = plan_DA.P_dso_grid(hour_idx_DA) + deltas.dso_grid(1);
    actual_operation.P_dso_seso_c(k)   = plan_DA.P_dso_charge(hour_idx_DA) + deltas.dso_seso_c(1);
    actual_operation.P_dso_seso_d(k)   = plan_DA.P_dso_discharge(hour_idx_DA) + deltas.dso_seso_d(1);
    
    actual_operation.P_mg1_net(k)      = plan_DA.P_mg1_net_exchange(hour_idx_DA) + deltas.mg1_net(1);
    actual_operation.P_mg1_seso_c(k)   = plan_DA.P_mg1_lease_c(hour_idx_DA) + deltas.mg1_seso_c(1);
    actual_operation.P_mg1_seso_d(k)   = plan_DA.P_mg1_lease_d(hour_idx_DA) + deltas.mg1_seso_d(1);

    actual_operation.P_mg2_net(k)      = plan_DA.P_mg2_net_exchange(hour_idx_DA) + deltas.mg2_net(1);
    actual_operation.P_mg2_seso_c(k)   = plan_DA.P_mg2_lease_c(hour_idx_DA) + deltas.mg2_seso_c(1);
    actual_operation.P_mg2_seso_d(k)   = plan_DA.P_mg2_lease_d(hour_idx_DA) + deltas.mg2_seso_d(1);

    actual_operation.P_mg3_net(k)      = plan_DA.P_mg3_net_exchange(hour_idx_DA) + deltas.mg3_net(1);
    actual_operation.P_mg3_seso_c(k)   = plan_DA.P_mg3_lease_c(hour_idx_DA) + deltas.mg3_seso_c(1);
    actual_operation.P_mg3_seso_d(k)   = plan_DA.P_mg3_lease_d(hour_idx_DA) + deltas.mg3_seso_d(1);
    
    % 记录详细的储能状态 (仅在成功时执行，失败时已在上面if语句中处理)
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
fprintf('日内滚动调度全部完成。\n\n');

%% 主函数：计及共享储能的配电网-多微网多时间尺度协同调度 (修正版)
clc; clear; close all;

fprintf('======================================================\n');
fprintf('      配电网-微网多时间尺度协同调度程序 V2.0 (集中式日内) \n');
fprintf('======================================================\n\n');

%% 1. 系统参数定义
T_DA = 1:24;
T_RT = (15/60):(15/60):24; % 修正时间轴以匹配96个点
H = 4;                     % 日内滚动优化窗口长度 (4个15分钟，即1小时)

%% 2. 数据准备
fprintf('(1/4) 正在准备日前及日内数据...\n');
dso_data = prepare_data_DSO(T_DA, T_RT);
mg_data = prepare_data_MGs(T_DA, T_RT);
fprintf('数据准备完毕。\n\n');

%% 3. 日前调度 (Day-Ahead Scheduling)
fprintf('(2/4) 正在执行日前经济调度...\n');
[plan_DA, admm_results_DA] = run_day_ahead_scheduling(dso_data.forecasts_DA, mg_data.forecasts_DA);
if isempty(plan_DA)
    error('日前调度失败，程序终止。');
end
fprintf('日前调度完成，已生成24小时运行计划。\n\n');

%% 4. 日内滚动调度 (Intra-day Rolling Scheduling)
fprintf('(3/4) 开始执行日内滚动调度 (总共 %d 个时间点)...\n', length(T_RT));

% 初始化用于存储日内实际运行结果的变量
actual_operation = struct();
actual_operation.P_dso_grid = zeros(1, 96);
actual_operation.P_dso_seso_c = zeros(1, 96);
actual_operation.P_dso_seso_d = zeros(1, 96);
actual_operation.P_mg1_net = zeros(1, 96);
actual_operation.P_mg2_net = zeros(1, 96);
actual_operation.P_mg3_net = zeros(1, 96);
actual_operation.P_mg1_seso_c = zeros(1, 96);
actual_operation.P_mg1_seso_d = zeros(1, 96);
actual_operation.P_mg2_seso_c = zeros(1, 96);
actual_operation.P_mg2_seso_d = zeros(1, 96);
actual_operation.P_mg3_seso_c = zeros(1, 96);
actual_operation.P_mg3_seso_d = zeros(1, 96);
actual_operation.E_seso = zeros(1, 96);
actual_operation.P_seso_ch = zeros(1, 96);
actual_operation.P_seso_dis = zeros(1, 96);
actual_operation.E_mg1 = zeros(1, 96);
actual_operation.P_mg1_batc = zeros(1, 96);
actual_operation.P_mg1_batd = zeros(1, 96);
actual_operation.E_mg2 = zeros(1, 96);
actual_operation.P_mg2_batc = zeros(1, 96);
actual_operation.P_mg2_batd = zeros(1, 96);
actual_operation.E_mg3 = zeros(1, 96);
actual_operation.P_mg3_batc = zeros(1, 96);
actual_operation.P_mg3_batd = zeros(1, 96);

% ======================== MODIFICATION START ========================
% 为 k=1 设置初始SOC (假设与日前计划一致或固定值)
initial_soc_mg1 = 800; % Or plan_DA.E_mg1(1);
initial_soc_mg2 = 800; % Or plan_DA.E_mg2(1);
initial_soc_mg3 = 800; % Or plan_DA.E_mg3(1);
% 注意：SESO的初始SOC由run_intraday_centralized内部从plan_DA读取，无需在此设置
% ======================== MODIFICATION END ========================


% 模拟日内滚动过程
for k = 1:96
    fprintf('--- 正在处理日内时间点 k = %d ---\n', k);
    
    % a. 确定当前滚动窗口
    window_idx = k:min(k + H - 1, 96);
    window_len = length(window_idx);

    % b. 获取当前窗口的实时预测数据
    current_forecasts_RT.dso = filter_data_by_window(dso_data.forecasts_RT, window_idx);
    current_forecasts_RT.mg  = filter_data_by_window(mg_data.forecasts_RT, window_idx);

    % ======================== MODIFICATION START ========================
    % c. 获取上一个时间点的实际SOC作为当前优化的初始SOC
    if k == 1
        current_initial_soc_mg1 = initial_soc_mg1;
        current_initial_soc_mg2 = initial_soc_mg2;
        current_initial_soc_mg3 = initial_soc_mg3;
    else
        % 从上一个时间点k-1的实际运行结果中获取SOC
        current_initial_soc_mg1 = actual_operation.E_mg1(k-1);
        current_initial_soc_mg2 = actual_operation.E_mg2(k-1);
        current_initial_soc_mg3 = actual_operation.E_mg3(k-1);
    end
    
    % d. 执行日内滚动优化 (传入初始SOC)
    [deltas, success] = run_intraday_centralized(plan_DA, current_forecasts_RT, k, window_len, ...
                                                 current_initial_soc_mg1, current_initial_soc_mg2, current_initial_soc_mg3);
    % ======================== MODIFICATION END ========================
                                                 
    % 获取当前时间点对应的日前小时索引
    hour_idx_DA = floor((k-1) / 4) + 1;
    
    if ~success || isempty(deltas)
        warning('在时间点 k = %d 的日内调度求解失败！将使用日前计划值 (调整量为0)。', k);
        deltas = get_zero_deltas(window_len);
        
        % 强制使用日前计划值填充储能状态 (注意：SOC需要特殊处理，因为是状态量)
        actual_operation.E_seso(k) = plan_DA.E_seso(hour_idx_DA); % 简化处理，实际应模拟演变
        actual_operation.P_seso_ch(k) = plan_DA.P_seso_ch(hour_idx_DA);
        actual_operation.P_seso_dis(k) = plan_DA.P_seso_dis(hour_idx_DA);
        
        % 对于MG储能，如果优化失败，SOC保持上一个时刻的值
        if k==1
             actual_operation.E_mg1(k) = initial_soc_mg1;
             actual_operation.E_mg2(k) = initial_soc_mg2;
             actual_operation.E_mg3(k) = initial_soc_mg3;
        else
             actual_operation.E_mg1(k) = actual_operation.E_mg1(k-1);
             actual_operation.E_mg2(k) = actual_operation.E_mg2(k-1);
             actual_operation.E_mg3(k) = actual_operation.E_mg3(k-1);
        end
        actual_operation.P_mg1_batc(k) = plan_DA.P_mg1_batc(hour_idx_DA); % 功率使用日前计划
        actual_operation.P_mg1_batd(k) = plan_DA.P_mg1_batd(hour_idx_DA);
        actual_operation.P_mg2_batc(k) = plan_DA.P_mg2_batc(hour_idx_DA);
        actual_operation.P_mg2_batd(k) = plan_DA.P_mg2_batd(hour_idx_DA);
        actual_operation.P_mg3_batc(k) = plan_DA.P_mg3_batc(hour_idx_DA);
        actual_operation.P_mg3_batd(k) = plan_DA.P_mg3_batd(hour_idx_DA);

    end

    % e. 记录第一个时间点的实际执行结果
    actual_operation.P_dso_grid(k)     = plan_DA.P_dso_grid(hour_idx_DA) + deltas.dso_grid(1);
    actual_operation.P_dso_seso_c(k)   = plan_DA.P_dso_charge(hour_idx_DA) + deltas.dso_seso_c(1);
    actual_operation.P_dso_seso_d(k)   = plan_DA.P_dso_discharge(hour_idx_DA) + deltas.dso_seso_d(1);
    
    actual_operation.P_mg1_net(k)      = plan_DA.P_mg1_net_exchange(hour_idx_DA) + deltas.mg1_net(1);
    actual_operation.P_mg1_seso_c(k)   = plan_DA.P_mg1_lease_c(hour_idx_DA) + deltas.mg1_seso_c(1);
    actual_operation.P_mg1_seso_d(k)   = plan_DA.P_mg1_lease_d(hour_idx_DA) + deltas.mg1_seso_d(1);

    actual_operation.P_mg2_net(k)      = plan_DA.P_mg2_net_exchange(hour_idx_DA) + deltas.mg2_net(1);
    actual_operation.P_mg2_seso_c(k)   = plan_DA.P_mg2_lease_c(hour_idx_DA) + deltas.mg2_seso_c(1);
    actual_operation.P_mg2_seso_d(k)   = plan_DA.P_mg2_lease_d(hour_idx_DA) + deltas.mg2_seso_d(1);

    actual_operation.P_mg3_net(k)      = plan_DA.P_mg3_net_exchange(hour_idx_DA) + deltas.mg3_net(1);
    actual_operation.P_mg3_seso_c(k)   = plan_DA.P_mg3_lease_c(hour_idx_DA) + deltas.mg3_seso_c(1);
    actual_operation.P_mg3_seso_d(k)   = plan_DA.P_mg3_lease_d(hour_idx_DA) + deltas.mg3_seso_d(1);
    
    % 记录详细的储能状态 (仅在成功时执行，失败时已在上面if语句中处理)
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
fprintf('日内滚动调度全部完成。\n\n');

% %% 5. 结果可视化 (V3: 全部改为单张 figure 输出)
% % ... (这是之前被注释掉的中文版可视化代码) ...
% fprintf('(4/4) 正在生成结果可视化图表...\n');
% time_axis_da = T_DA;
% time_axis_rt = T_RT;
% 
% % 图 1: 日前优化ADMM收敛曲线 - 目标函数
% figure('Name', '日前优化ADMM目标函数收敛性');
% plot(1:length(admm_results_DA.Obj_DSO), admm_results_DA.Obj_DSO, 'r-o', 'LineWidth', 1.5); hold on;
% plot(1:length(admm_results_DA.Obj_MG1), admm_results_DA.Obj_MG1, 'g-s', 'LineWidth', 1.5);
% plot(1:length(admm_results_DA.Obj_MG2), admm_results_DA.Obj_MG2, 'b-^', 'LineWidth', 1.5);
% plot(1:length(admm_results_DA.Obj_MG3), admm_results_DA.Obj_MG3, 'm-d', 'LineWidth', 1.5);
% plot(1:length(admm_results_DA.Obj_SESO), admm_results_DA.Obj_SESO, 'k-*', 'LineWidth', 1.5);
% title('日前优化各主体目标函数收敛过程');
% xlabel('迭代次数'); ylabel('目标函数值');
% legend('DSO', 'MG1', 'MG2', 'MG3', 'SESO', 'Location', 'best');
% % grid on; % <--- 修改
% 
% % 图 2: 日前优化ADMM收敛曲线 - 残差
% figure('Name', '日前优化ADMM残差收敛性');
% semilogy(1:length(admm_results_DA.Residual), admm_results_DA.Residual, 'b-x', 'LineWidth', 1.5);
% title('日前优化ADMM算法残差收敛过程');
% xlabel('迭代次数'); ylabel('残差 (对数尺度)');
% % grid on; % <--- 修改
% 
% % --- 新能源出力对比 ---
% dso_re_da = sum(dso_data.forecasts_DA.P_pv, 1) + sum(dso_data.forecasts_DA.P_wt, 1);
% dso_re_rt = sum(dso_data.forecasts_RT.P_pv, 1) + sum(dso_data.forecasts_RT.P_wt, 1);
% figure('Name', 'DSO新能源出力对比');
% stairs(time_axis_da, dso_re_da, 'r--', 'LineWidth', 1.5); hold on;
% plot(time_axis_rt, dso_re_rt, 'b-', 'LineWidth', 1.5);
% title('DSO 区域总新能源出力'); xlabel('时间 (h)'); ylabel('功率 (kW)');
% legend('日前预测', '日内实际', 'Location', 'best'); % grid on; % <--- 修改
% xlim([0, 24]);
% 
% figure('Name', 'MG1新能源出力对比');
% stairs(time_axis_da, mg_data.forecasts_DA.Predict_wt_mg1, 'r--', 'LineWidth', 1.5); hold on;
% plot(time_axis_rt, mg_data.forecasts_RT.Predict_wt_mg1, 'b-', 'LineWidth', 1.5);
% title('微网1 风电出力'); xlabel('时间 (h)'); ylabel('功率 (kW)');
% legend('日前预测', '日内实际', 'Location', 'best'); % grid on; % <--- 修改
% xlim([0, 24]);
% 
% figure('Name', 'MG2新能源出力对比');
% stairs(time_axis_da, mg_data.forecasts_DA.Predict_pv_mg2, 'r--', 'LineWidth', 1.5); hold on;
% plot(time_axis_rt, mg_data.forecasts_RT.Predict_pv_mg2, 'b-', 'LineWidth', 1.5);
% title('微网2 光伏出力'); xlabel('时间 (h)'); ylabel('功率 (kW)');
% legend('日前预测', '日内实际', 'Location', 'best'); % grid on; % <--- 修改
% xlim([0, 24]);
% 
% figure('Name', 'MG3新能源出力对比');
% stairs(time_axis_da, mg_data.forecasts_DA.Predict_pv_mg3, 'r--', 'LineWidth', 1.5); hold on;
% plot(time_axis_rt, mg_data.forecasts_RT.Predict_pv_mg3, 'b-', 'LineWidth', 1.5);
% title('微网3 光伏出力'); xlabel('时间 (h)'); ylabel('功率 (kW)');
% legend('日前预测', '日内实际', 'Location', 'best'); % grid on; % <--- 修改
% xlim([0, 24]);
% 
% % --- 负荷对比 ---
% figure('Name', 'DSO负荷对比');
% stairs(time_axis_da, sum(dso_data.forecasts_DA.PLOAD, 1), 'r--', 'LineWidth', 1.5); hold on;
% plot(time_axis_rt, sum(dso_data.forecasts_RT.PLOAD, 1), 'b-', 'LineWidth', 1.5);
% title('DSO 区域总负荷'); xlabel('时间 (h)'); ylabel('功率 (kW)');
% legend('日前预测', '日内实际', 'Location', 'best'); % grid on; % <--- 修改
% xlim([0, 24]);
% 
% figure('Name', 'MG1负荷对比');
% stairs(time_axis_da, mg_data.forecasts_DA.L_e0_mg1, 'r--', 'LineWidth', 1.5); hold on;
% plot(time_axis_rt, mg_data.forecasts_RT.L_e0_mg1, 'b-', 'LineWidth', 1.5);
% title('微网1 电负荷'); xlabel('时间 (h)'); ylabel('功率 (kW)');
% legend('日前预测', '日内实际', 'Location', 'best'); % grid on; % <--- 修改
% xlim([0, 24]);
% 
% figure('Name', 'MG2负荷对比');
% stairs(time_axis_da, mg_data.forecasts_DA.L_e0_mg2, 'r--', 'LineWidth', 1.5); hold on;
% plot(time_axis_rt, mg_data.forecasts_RT.L_e0_mg2, 'b-', 'LineWidth', 1.5);
% title('微网2 电负荷'); xlabel('时间 (h)'); ylabel('功率 (kW)');
% legend('日前预测', '日内实际', 'Location', 'best'); % grid on; % <--- 修改
% xlim([0, 24]);
% 
% figure('Name', 'MG3负荷对比');
% stairs(time_axis_da, mg_data.forecasts_DA.L_e0_mg3, 'r--', 'LineWidth', 1.5); hold on;
% plot(time_axis_rt, mg_data.forecasts_RT.L_e0_mg3, 'b-', 'LineWidth', 1.5);
% title('微网3 电负荷'); xlabel('时间 (h)'); ylabel('功率 (kW)');
% legend('日前预测', '日内实际', 'Location', 'best'); % grid on; % <--- 修改
% xlim([0, 24]);
% 
% % --- 储能租赁情况对比 (充电) ---
% figure('Name', 'DSO向SESO充电');
% stairs(time_axis_da, plan_DA.P_dso_charge, 'r--', 'LineWidth', 1.5); hold on;
% plot(time_axis_rt, actual_operation.P_dso_seso_c, 'b-', 'LineWidth', 1.5);
% title('DSO向SESO充电'); xlabel('时间 (h)'); ylabel('功率 (kW)');
% legend('日前计划', '日内实际', 'Location', 'best'); % grid on; % <--- 修改
% xlim([0, 24]);
% 
% figure('Name', 'MG1向SESO提供充电服务');
% stairs(time_axis_da, plan_DA.P_mg1_lease_c, 'r--', 'LineWidth', 1.5); hold on;
% plot(time_axis_rt, actual_operation.P_mg1_seso_c, 'b-', 'LineWidth', 1.5);
% title('MG1向SESO提供充电服务'); xlabel('时间 (h)'); ylabel('功率 (kW)');
% legend('日前计划', '日内实际', 'Location', 'best'); % <--- 修改
% % grid on; % <--- 修改
% xlim([0, 24]);
% 
% figure('Name', 'MG2向SESO提供充电服务');
% stairs(time_axis_da, plan_DA.P_mg2_lease_c, 'r--', 'LineWidth', 1.5); hold on;
% plot(time_axis_rt, actual_operation.P_mg2_seso_c, 'b-', 'LineWidth', 1.5);
% title('MG2向SESO提供充电服务'); xlabel('时间 (h)'); ylabel('功率 (kW)');
% legend('日前计划', '日内实际', 'Location', 'best'); % <--- 修改
% % grid on; % <--- 修改
% xlim([0, 24]);
% 
% figure('Name', 'MG3向SESO提供充电服务');
% stairs(time_axis_da, plan_DA.P_mg3_lease_c, 'r--', 'LineWidth', 1.5); hold on;
% plot(time_axis_rt, actual_operation.P_mg3_seso_c, 'b-', 'LineWidth', 1.5);
% title('MG3向SESO提供充电服务'); xlabel('时间 (h)'); ylabel('功率 (kW)');
% legend('日前计划', '日内实际', 'Location', 'best'); % <--- 修改
% % grid on; % <--- 修改
% xlim([0, 24]);
% 
% % --- 储能租赁情况对比 (放电) ---
% figure('Name', 'DSO向SESO放电');
% stairs(time_axis_da, plan_DA.P_dso_discharge, 'r--', 'LineWidth', 1.5); hold on;
% plot(time_axis_rt, actual_operation.P_dso_seso_d, 'b-', 'LineWidth', 1.5);
% title('DSO向SESO放电'); xlabel('时间 (h)'); ylabel('功率 (kW)');
% legend('日前计划', '日内实际', 'Location', 'best'); % grid on; % <--- 修改
% xlim([0, 24]);
% 
% figure('Name', 'MG1向SESO提供放电服务');
% stairs(time_axis_da, plan_DA.P_mg1_lease_d, 'r--', 'LineWidth', 1.5); hold on;
% plot(time_axis_rt, actual_operation.P_mg1_seso_d, 'b-', 'LineWidth', 1.5);
% title('MG1向SESO提供放电服务'); xlabel('时间 (h)'); ylabel('功率 (kW)');
% legend('日前计划', '日内实际', 'Location', 'best'); % <--- 修改
% % grid on; % <--- 修改
% xlim([0, 24]);
% 
% figure('Name', 'MG2向SESO提供放电服务');
% stairs(time_axis_da, plan_DA.P_mg2_lease_d, 'r--', 'LineWidth', 1.5); hold on;
% plot(time_axis_rt, actual_operation.P_mg2_seso_d, 'b-', 'LineWidth', 1.5);
% title('MG2向SESO提供放电服务'); xlabel('时间 (h)'); ylabel('功率 (kW)');
% legend('日前计划', '日内实际', 'Location', 'best'); % <--- 修改
% % grid on; % <--- 修改
% xlim([0, 24]);
% 
% figure('Name', 'MG3向SESO提供放电服务');
% stairs(time_axis_da, plan_DA.P_mg3_lease_d, 'r--', 'LineWidth', 1.5); hold on;
% plot(time_axis_rt, actual_operation.P_mg3_seso_d, 'b-', 'LineWidth', 1.5);
% title('MG3向SESO提供放电服务'); xlabel('时间 (h)'); ylabel('功率 (kW)');
% legend('日前计划', '日内实际', 'Location', 'best'); % <--- 修改
% % grid on; % <--- 修改
% xlim([0, 24]);
% 
% % --- SESO 自有储能状态对比 ---
% figure('Name', 'SESO 自有储能SOC对比');
% stairs(time_axis_da, plan_DA.E_seso(1:24), 'r--', 'LineWidth', 1.5); hold on;
% plot(time_axis_rt, actual_operation.E_seso, 'b-', 'LineWidth', 1.5);
% title('SESO 自有储能SOC (日前 vs 日内)');
% xlabel('时间 (h)'); ylabel('电量 (kWh)');
% legend('日前计划', '日内实际', 'Location', 'best');
% % grid on;
% xlim([0, 24]);
% 
% figure('Name', 'SESO 自有储能充放电功率对比');
% stairs(time_axis_da, plan_DA.P_seso_ch, 'r--', 'LineWidth', 1.5); hold on;
% plot(time_axis_rt, actual_operation.P_seso_ch, 'b-', 'LineWidth', 1.5);
% stairs(time_axis_da, -plan_DA.P_seso_dis, 'm--', 'LineWidth', 1.5);
% plot(time_axis_rt, -actual_operation.P_seso_dis, 'c-', 'LineWidth', 1.5);
% yline(0, 'k:', 'LineWidth', 1);
% title('SESO 自有储能充放电功率 (日前 vs 日内)');
% xlabel('时间 (h)'); ylabel('功率 (kW)');
% legend('日前充电', '日内充电', '日前放电', '日内放电', 'Location', 'best');
% % grid on;
% xlim([0, 24]);
% 
% % --- MG1 自有储能状态对比 ---
% figure('Name', 'MG1 自有储能SOC对比');
% stairs(time_axis_da, plan_DA.E_mg1, 'r--', 'LineWidth', 1.5); hold on;
% plot(time_axis_rt, actual_operation.E_mg1, 'b-', 'LineWidth', 1.5);
% title('MG1 自有储能SOC (日前 vs 日内)');
% xlabel('时间 (h)'); ylabel('电量 (kWh)');
% legend('日前计划', '日内实际', 'Location', 'best');
% % grid on;
% xlim([0, 24]);
% 
% figure('Name', 'MG1 自有储能充放电功率对比');
% stairs(time_axis_da, plan_DA.P_mg1_batc, 'r--', 'LineWidth', 1.5); hold on;
% plot(time_axis_rt, actual_operation.P_mg1_batc, 'b-', 'LineWidth', 1.5);
% stairs(time_axis_da, -plan_DA.P_mg1_batd, 'm--', 'LineWidth', 1.5);
% plot(time_axis_rt, -actual_operation.P_mg1_batd, 'c-', 'LineWidth', 1.5);
% yline(0, 'k:', 'LineWidth', 1);
% title('MG1 自有储能充放电功率 (日前 vs 日内)');
% xlabel('时间 (h)'); ylabel('功率 (kW)');
% legend('日前充电', '日内充电', '日前放电', '日内放电', 'Location', 'best');
% % grid on;
% xlim([0, 24]);
% 
% % --- MG2 自有储能状态对比 ---
% figure('Name', 'MG2 自有储能SOC对比');
% stairs(time_axis_da, plan_DA.E_mg2, 'r--', 'LineWidth', 1.5); hold on;
% plot(time_axis_rt, actual_operation.E_mg2, 'b-', 'LineWidth', 1.5);
% title('MG2 自有储能SOC (日前 vs 日内)');
% xlabel('时间 (h)'); ylabel('电量 (kWh)');
% legend('日前计划', '日内实际', 'Location', 'best');
% % grid on;
% xlim([0, 24]);
% 
% figure('Name', 'MG2 自有储能充放电功率对比');
% stairs(time_axis_da, plan_DA.P_mg2_batc, 'r--', 'LineWidth', 1.5); hold on;
% plot(time_axis_rt, actual_operation.P_mg2_batc, 'b-', 'LineWidth', 1.5);
% stairs(time_axis_da, -plan_DA.P_mg2_batd, 'm--', 'LineWidth', 1.5);
% plot(time_axis_rt, -actual_operation.P_mg2_batd, 'c-', 'LineWidth', 1.5);
% yline(0, 'k:', 'LineWidth', 1);
% title('MG2 自有储能充放电功率 (日前 vs 日内)');
% xlabel('时间 (h)'); ylabel('功率 (kW)');
% legend('日前充电', '日内充电', '日前放电', '日内放电', 'Location', 'best');
% % grid on;
% xlim([0, 24]);
% 
% % --- MG3 自有储能状态对比 ---
% figure('Name', 'MG3 自有储能SOC对比');
% stairs(time_axis_da, plan_DA.E_mg3, 'r--', 'LineWidth', 1.5); hold on;
% plot(time_axis_rt, actual_operation.E_mg3, 'b-', 'LineWidth', 1.5);
% title('MG3 自有储能SOC (日前 vs 日内)');
% xlabel('时间 (h)'); ylabel('电量 (kWh)');
% legend('日前计划', '日内实际', 'Location', 'best');
% % grid on;
% xlim([0, 24]);
% 
% figure('Name', 'MG3 自有储能充放电功率对比');
% stairs(time_axis_da, plan_DA.P_mg3_batc, 'r--', 'LineWidth', 1.5); hold on;
% plot(time_axis_rt, actual_operation.P_mg3_batc, 'b-', 'LineWidth', 1.5);
% stairs(time_axis_da, -plan_DA.P_mg3_batd, 'm--', 'LineWidth', 1.5);
% plot(time_axis_rt, -actual_operation.P_mg3_batd, 'c-', 'LineWidth', 1.5);
% yline(0, 'k:', 'LineWidth', 1);
% title('MG3 自有储能充放电功率 (日前 vs 日内)');
% xlabel('时间 (h)'); ylabel('功率 (kW)');
% legend('日前充电', '日内充电', '日前放电', '日内放电', 'Location', 'best');
% % grid on;
% xlim([0, 24]);
% 
% 
% fprintf('程序运行结束。\n');

%% 5. 结果可视化 (V4: 英文版, 字体放大, Y轴扩展, 图例固定)
% Visualization (V4: English Version, Fonts Enlarged, Y-Axis Expanded, Legend Fixed)
fprintf('(4/4) Generating result visualization charts...\n');
time_axis_da = T_DA;
time_axis_rt = T_RT;
base_font_size = 14;  % MODIFIED: Increased font size
title_font_size = 16; % MODIFIED: Increased font size

% --- 图 1: 日前优化ADMM收敛曲线 - 目标函数 ---
% Figure 1: Day-Ahead ADMM Convergence - Objective Function
figure('Name', 'DA ADMM Objective Convergence');
plot(1:length(admm_results_DA.Obj_DSO), admm_results_DA.Obj_DSO, 'r-o', 'LineWidth', 1.5); hold on;
plot(1:length(admm_results_DA.Obj_MG1), admm_results_DA.Obj_MG1, 'g-s', 'LineWidth', 1.5);
plot(1:length(admm_results_DA.Obj_MG2), admm_results_DA.Obj_MG2, 'b-^', 'LineWidth', 1.5);
plot(1:length(admm_results_DA.Obj_MG3), admm_results_DA.Obj_MG3, 'm-d', 'LineWidth', 1.5);
plot(1:length(admm_results_DA.Obj_SESO), admm_results_DA.Obj_SESO, 'k-*', 'LineWidth', 1.5);
title('Day-Ahead ADMM Objective Function Convergence', 'FontSize', title_font_size);
xlabel('Iterations', 'FontSize', base_font_size); 
ylabel('Objective Value', 'FontSize', base_font_size);
lgd = legend('DSO', 'MG1', 'MG2', 'MG3', 'SESO', 'Location', 'northeast'); % MODIFIED: Location
set(lgd, 'FontSize', base_font_size);
set(gca, 'FontSize', base_font_size);
% MODIFIED: Add code to expand Y-axis
ax = gca; ylim_vals = ax.YLim;
if ylim_vals(2) < 0.1, ylim_vals(2) = 1; end;
if ylim_vals(1) > -0.1, ylim_vals(1) = -1; end;
ax.YLim = [ylim_vals(1) * 1.15, ylim_vals(2) * 1.15];

% --- 图 2: 日前优化ADMM收敛曲线 - 残差 ---
% Figure 2: Day-Ahead ADMM Convergence - Residual
figure('Name', 'DA ADMM Residual Convergence');
semilogy(1:length(admm_results_DA.Residual), admm_results_DA.Residual, 'b-x', 'LineWidth', 1.5);
title('Day-Ahead ADMM Residual Convergence', 'FontSize', title_font_size);
xlabel('Iterations', 'FontSize', base_font_size); 
ylabel('Residual (Log Scale)', 'FontSize', base_font_size);
set(gca, 'FontSize', base_font_size);
% MODIFIED: Add code to expand Y-axis (Top only for log)
ax = gca; ylim_vals = ax.YLim;
if ylim_vals(2) < 1, ylim_vals(2) = 1; end;
ax.YLim = [ylim_vals(1), ylim_vals(2) * 1.15];

% --- 新能源出力对比 ---
% --- RE Output Comparison ---
dso_re_da = sum(dso_data.forecasts_DA.P_pv, 1) + sum(dso_data.forecasts_DA.P_wt, 1);
dso_re_rt = sum(dso_data.forecasts_RT.P_pv, 1) + sum(dso_data.forecasts_RT.P_wt, 1);
figure('Name', 'DSO RE Output Comparison');
stairs(time_axis_da, dso_re_da, 'r--', 'LineWidth', 1.5); hold on;
plot(time_axis_rt, dso_re_rt, 'b-', 'LineWidth', 1.5);
title('DSO Total RE Output', 'FontSize', title_font_size); 
xlabel('Time (h)', 'FontSize', base_font_size); 
ylabel('Power (MW)', 'FontSize', base_font_size); % Note: DSO data is often in MW
lgd = legend('Day-Ahead Forecast', 'Intra-day Actual', 'Location', 'northeast'); % MODIFIED: Location
set(lgd, 'FontSize', base_font_size);
set(gca, 'FontSize', base_font_size);
xlim([0, 24]);
% MODIFIED: Add code to expand Y-axis (Top only)
ax = gca; ylim_vals = ax.YLim;
if ylim_vals(2) < 0.1, ylim_vals(2) = 1; end;
if ylim_vals(1) > -0.1 && ylim_vals(1) < 0.1, ax.YLim = [0, ylim_vals(2) * 1.15]; else, ax.YLim = [ylim_vals(1), ylim_vals(2) * 1.15]; end

figure('Name', 'MG1 RE Output Comparison');
stairs(time_axis_da, mg_data.forecasts_DA.Predict_wt_mg1, 'r--', 'LineWidth', 1.5); hold on;
plot(time_axis_rt, mg_data.forecasts_RT.Predict_wt_mg1, 'b-', 'LineWidth', 1.5);
title('Microgrid 1 Wind Power Output', 'FontSize', title_font_size); 
xlabel('Time (h)', 'FontSize', base_font_size); 
ylabel('Power (kW)', 'FontSize', base_font_size);
lgd = legend('Day-Ahead Forecast', 'Intra-day Actual', 'Location', 'northeast'); % MODIFIED: Location
set(lgd, 'FontSize', base_font_size);
set(gca, 'FontSize', base_font_size);
xlim([0, 24]);
% MODIFIED: Add code to expand Y-axis (Top only)
ax = gca; ylim_vals = ax.YLim;
if ylim_vals(2) < 0.1, ylim_vals(2) = 1; end;
if ylim_vals(1) > -0.1 && ylim_vals(1) < 0.1, ax.YLim = [0, ylim_vals(2) * 1.15]; else, ax.YLim = [ylim_vals(1), ylim_vals(2) * 1.15]; end

figure('Name', 'MG2 RE Output Comparison');
stairs(time_axis_da, mg_data.forecasts_DA.Predict_pv_mg2, 'r--', 'LineWidth', 1.5); hold on;
plot(time_axis_rt, mg_data.forecasts_RT.Predict_pv_mg2, 'b-', 'LineWidth', 1.5);
title('Microgrid 2 PV Output', 'FontSize', title_font_size); 
xlabel('Time (h)', 'FontSize', base_font_size); 
ylabel('Power (kW)', 'FontSize', base_font_size);
lgd = legend('Day-Ahead Forecast', 'Intra-day Actual', 'Location', 'northeast'); % MODIFIED: Location
set(lgd, 'FontSize', base_font_size);
set(gca, 'FontSize', base_font_size);
xlim([0, 24]);
% MODIFIED: Add code to expand Y-axis (Top only)
ax = gca; ylim_vals = ax.YLim;
if ylim_vals(2) < 0.1, ylim_vals(2) = 1; end;
if ylim_vals(1) > -0.1 && ylim_vals(1) < 0.1, ax.YLim = [0, ylim_vals(2) * 1.15]; else, ax.YLim = [ylim_vals(1), ylim_vals(2) * 1.15]; end

figure('Name', 'MG3 RE Output Comparison');
stairs(time_axis_da, mg_data.forecasts_DA.Predict_pv_mg3, 'r--', 'LineWidth', 1.5); hold on;
plot(time_axis_rt, mg_data.forecasts_RT.Predict_pv_mg3, 'b-', 'LineWidth', 1.5);
title('Microgrid 3 PV Output', 'FontSize', title_font_size); 
xlabel('Time (h)', 'FontSize', base_font_size); 
ylabel('Power (kW)', 'FontSize', base_font_size);
lgd = legend('Day-Ahead Forecast', 'Intra-day Actual', 'Location', 'northeast'); % MODIFIED: Location
set(lgd, 'FontSize', base_font_size);
set(gca, 'FontSize', base_font_size);
xlim([0, 24]);
% MODIFIED: Add code to expand Y-axis (Top only)
ax = gca; ylim_vals = ax.YLim;
if ylim_vals(2) < 0.1, ylim_vals(2) = 1; end;
if ylim_vals(1) > -0.1 && ylim_vals(1) < 0.1, ax.YLim = [0, ylim_vals(2) * 1.15]; else, ax.YLim = [ylim_vals(1), ylim_vals(2) * 1.15]; end

% --- 负荷对比 ---
% --- Load Comparison ---
figure('Name', 'DSO Load Comparison');
stairs(time_axis_da, sum(dso_data.forecasts_DA.PLOAD, 1), 'r--', 'LineWidth', 1.5); hold on;
plot(time_axis_rt, sum(dso_data.forecasts_RT.PLOAD, 1), 'b-', 'LineWidth', 1.5);
title('DSO Total Load', 'FontSize', title_font_size); 
xlabel('Time (h)', 'FontSize', base_font_size); 
ylabel('Power (kW)', 'FontSize', base_font_size);
lgd = legend('Day-Ahead Forecast', 'Intra-day Actual', 'Location', 'northeast'); % MODIFIED: Location
set(lgd, 'FontSize', base_font_size);
set(gca, 'FontSize', base_font_size);
xlim([0, 24]);
% MODIFIED: Add code to expand Y-axis (Top only)
ax = gca; ylim_vals = ax.YLim;
if ylim_vals(2) < 0.1, ylim_vals(2) = 1; end;
if ylim_vals(1) > -0.1 && ylim_vals(1) < 0.1, ax.YLim = [0, ylim_vals(2) * 1.15]; else, ax.YLim = [ylim_vals(1), ylim_vals(2) * 1.15]; end

figure('Name', 'MG1 Load Comparison');
stairs(time_axis_da, mg_data.forecasts_DA.L_e0_mg1, 'r--', 'LineWidth', 1.5); hold on;
plot(time_axis_rt, mg_data.forecasts_RT.L_e0_mg1, 'b-', 'LineWidth', 1.5);
title('Microgrid 1 Electrical Load', 'FontSize', title_font_size); 
xlabel('Time (h)', 'FontSize', base_font_size); 
ylabel('Power (kW)', 'FontSize', base_font_size);
lgd = legend('Day-Ahead Forecast', 'Intra-day Actual', 'Location', 'northeast'); % MODIFIED: Location
set(lgd, 'FontSize', base_font_size);
set(gca, 'FontSize', base_font_size);
xlim([0, 24]);
% MODIFIED: Add code to expand Y-axis (Top only)
ax = gca; ylim_vals = ax.YLim;
if ylim_vals(2) < 0.1, ylim_vals(2) = 1; end;
if ylim_vals(1) > -0.1 && ylim_vals(1) < 0.1, ax.YLim = [0, ylim_vals(2) * 1.15]; else, ax.YLim = [ylim_vals(1), ylim_vals(2) * 1.15]; end

figure('Name', 'MG2 Load Comparison');
stairs(time_axis_da, mg_data.forecasts_DA.L_e0_mg2, 'r--', 'LineWidth', 1.5); hold on;
plot(time_axis_rt, mg_data.forecasts_RT.L_e0_mg2, 'b-', 'LineWidth', 1.5);
title('Microgrid 2 Electrical Load', 'FontSize', title_font_size); 
xlabel('Time (h)', 'FontSize', base_font_size); 
ylabel('Power (kW)', 'FontSize', base_font_size);
lgd = legend('Day-Ahead Forecast', 'Intra-day Actual', 'Location', 'northeast'); % MODIFIED: Location
set(lgd, 'FontSize', base_font_size);
set(gca, 'FontSize', base_font_size);
xlim([0, 24]);
% MODIFIED: Add code to expand Y-axis (Top only)
ax = gca; ylim_vals = ax.YLim;
if ylim_vals(2) < 0.1, ylim_vals(2) = 1; end;
if ylim_vals(1) > -0.1 && ylim_vals(1) < 0.1, ax.YLim = [0, ylim_vals(2) * 1.15]; else, ax.YLim = [ylim_vals(1), ylim_vals(2) * 1.15]; end

figure('Name', 'MG3 Load Comparison');
stairs(time_axis_da, mg_data.forecasts_DA.L_e0_mg3, 'r--', 'LineWidth', 1.5); hold on;
plot(time_axis_rt, mg_data.forecasts_RT.L_e0_mg3, 'b-', 'LineWidth', 1.5);
title('Microgrid 3 Electrical Load', 'FontSize', title_font_size); 
xlabel('Time (h)', 'FontSize', base_font_size); 
ylabel('Power (kW)', 'FontSize', base_font_size);
lgd = legend('Day-Ahead Forecast', 'Intra-day Actual', 'Location', 'northeast'); % MODIFIED: Location
set(lgd, 'FontSize', base_font_size);
set(gca, 'FontSize', base_font_size);
xlim([0, 24]);
% MODIFIED: Add code to expand Y-axis (Top only)
ax = gca; ylim_vals = ax.YLim;
if ylim_vals(2) < 0.1, ylim_vals(2) = 1; end;
if ylim_vals(1) > -0.1 && ylim_vals(1) < 0.1, ax.YLim = [0, ylim_vals(2) * 1.15]; else, ax.YLim = [ylim_vals(1), ylim_vals(2) * 1.15]; end

% --- 储能租赁情况对比 (充电) ---
% --- Storage Leasing Comparison (Charging) ---
figure('Name', 'DSO to SESO Charging');
stairs(time_axis_da, plan_DA.P_dso_charge, 'r--', 'LineWidth', 1.5); hold on;
plot(time_axis_rt, actual_operation.P_dso_seso_c, 'b-', 'LineWidth', 1.5);
title('DSO Charging SESO', 'FontSize', title_font_size); 
xlabel('Time (h)', 'FontSize', base_font_size); 
ylabel('Power (kW)', 'FontSize', base_font_size);
lgd = legend('Day-Ahead Plan', 'Intra-day Actual', 'Location', 'northeast'); % MODIFIED: Location
set(lgd, 'FontSize', base_font_size);
set(gca, 'FontSize', base_font_size);
xlim([0, 24]);
% MODIFIED: Add code to expand Y-axis (Top only)
ax = gca; ylim_vals = ax.YLim;
if ylim_vals(2) < 0.1, ylim_vals(2) = 1; end;
if ylim_vals(1) > -0.1 && ylim_vals(1) < 0.1, ax.YLim = [0, ylim_vals(2) * 1.15]; else, ax.YLim = [ylim_vals(1), ylim_vals(2) * 1.15]; end

figure('Name', 'MG1 Charging Service to SESO');
stairs(time_axis_da, plan_DA.P_mg1_lease_c, 'r--', 'LineWidth', 1.5); hold on;
plot(time_axis_rt, actual_operation.P_mg1_seso_c, 'b-', 'LineWidth', 1.5);
title('MG1 Charging Service to SESO', 'FontSize', title_font_size); 
xlabel('Time (h)', 'FontSize', base_font_size); 
ylabel('Power (kW)', 'FontSize', base_font_size);
lgd = legend('Day-Ahead Plan', 'Intra-day Actual', 'Location', 'northeast'); % MODIFIED: Location
set(lgd, 'FontSize', base_font_size);
set(gca, 'FontSize', base_font_size);
xlim([0, 24]);
% MODIFIED: Add code to expand Y-axis (Top only)
ax = gca; ylim_vals = ax.YLim;
if ylim_vals(2) < 0.1, ylim_vals(2) = 1; end;
if ylim_vals(1) > -0.1 && ylim_vals(1) < 0.1, ax.YLim = [0, ylim_vals(2) * 1.15]; else, ax.YLim = [ylim_vals(1), ylim_vals(2) * 1.15]; end

figure('Name', 'MG2 Charging Service to SESO');
stairs(time_axis_da, plan_DA.P_mg2_lease_c, 'r--', 'LineWidth', 1.5); hold on;
plot(time_axis_rt, actual_operation.P_mg2_seso_c, 'b-', 'LineWidth', 1.5);
title('MG2 Charging Service to SESO', 'FontSize', title_font_size); 
xlabel('Time (h)', 'FontSize', base_font_size); 
ylabel('Power (kW)', 'FontSize', base_font_size);
lgd = legend('Day-Ahead Plan', 'Intra-day Actual', 'Location', 'northeast'); % MODIFIED: Location
set(lgd, 'FontSize', base_font_size);
set(gca, 'FontSize', base_font_size);
xlim([0, 24]);
% MODIFIED: Add code to expand Y-axis (Top only)
ax = gca; ylim_vals = ax.YLim;
if ylim_vals(2) < 0.1, ylim_vals(2) = 1; end;
if ylim_vals(1) > -0.1 && ylim_vals(1) < 0.1, ax.YLim = [0, ylim_vals(2) * 1.15]; else, ax.YLim = [ylim_vals(1), ylim_vals(2) * 1.15]; end

figure('Name', 'MG3 Charging Service to SESO');
stairs(time_axis_da, plan_DA.P_mg3_lease_c, 'r--', 'LineWidth', 1.5); hold on;
plot(time_axis_rt, actual_operation.P_mg3_seso_c, 'b-', 'LineWidth', 1.5);
title('MG3 Charging Service to SESO', 'FontSize', title_font_size); 
xlabel('Time (h)', 'FontSize', base_font_size); 
ylabel('Power (kW)', 'FontSize', base_font_size);
lgd = legend('Day-Ahead Plan', 'Intra-day Actual', 'Location', 'northeast'); % MODIFIED: Location
set(lgd, 'FontSize', base_font_size);
set(gca, 'FontSize', base_font_size);
xlim([0, 24]);
% MODIFIED: Add code to expand Y-axis (Top only)
ax = gca; ylim_vals = ax.YLim;
if ylim_vals(2) < 0.1, ylim_vals(2) = 1; end;
if ylim_vals(1) > -0.1 && ylim_vals(1) < 0.1, ax.YLim = [0, ylim_vals(2) * 1.15]; else, ax.YLim = [ylim_vals(1), ylim_vals(2) * 1.15]; end

% --- 储能租赁情况对比 (放电) ---
% --- Storage Leasing Comparison (Discharging) ---
figure('Name', 'DSO from SESO Discharging');
stairs(time_axis_da, plan_DA.P_dso_discharge, 'r--', 'LineWidth', 1.5); hold on;
plot(time_axis_rt, actual_operation.P_dso_seso_d, 'b-', 'LineWidth', 1.5);
title('DSO Discharging from SESO', 'FontSize', title_font_size); 
xlabel('Time (h)', 'FontSize', base_font_size); 
ylabel('Power (kW)', 'FontSize', base_font_size);
lgd = legend('Day-Ahead Plan', 'Intra-day Actual', 'Location', 'northeast'); % MODIFIED: Location
set(lgd, 'FontSize', base_font_size);
set(gca, 'FontSize', base_font_size);
xlim([0, 24]);
% MODIFIED: Add code to expand Y-axis (Top only)
ax = gca; ylim_vals = ax.YLim;
if ylim_vals(2) < 0.1, ylim_vals(2) = 1; end;
if ylim_vals(1) > -0.1 && ylim_vals(1) < 0.1, ax.YLim = [0, ylim_vals(2) * 1.15]; else, ax.YLim = [ylim_vals(1), ylim_vals(2) * 1.15]; end

figure('Name', 'MG1 Discharging Service to SESO');
stairs(time_axis_da, plan_DA.P_mg1_lease_d, 'r--', 'LineWidth', 1.5); hold on;
plot(time_axis_rt, actual_operation.P_mg1_seso_d, 'b-', 'LineWidth', 1.5);
title('MG1 Discharging Service to SESO', 'FontSize', title_font_size); 
xlabel('Time (h)', 'FontSize', base_font_size); 
ylabel('Power (kW)', 'FontSize', base_font_size);
lgd = legend('Day-Ahead Plan', 'Intra-day Actual', 'Location', 'northeast'); % MODIFIED: Location
set(lgd, 'FontSize', base_font_size);
set(gca, 'FontSize', base_font_size);
xlim([0, 24]);
% MODIFIED: Add code to expand Y-axis (Top only)
ax = gca; ylim_vals = ax.YLim;
if ylim_vals(2) < 0.1, ylim_vals(2) = 1; end;
if ylim_vals(1) > -0.1 && ylim_vals(1) < 0.1, ax.YLim = [0, ylim_vals(2) * 1.15]; else, ax.YLim = [ylim_vals(1), ylim_vals(2) * 1.15]; end

figure('Name', 'MG2 Discharging Service to SESO');
stairs(time_axis_da, plan_DA.P_mg2_lease_d, 'r--', 'LineWidth', 1.5); hold on;
plot(time_axis_rt, actual_operation.P_mg2_seso_d, 'b-', 'LineWidth', 1.5);
title('MG2 Discharging Service to SESO', 'FontSize', title_font_size); 
xlabel('Time (h)', 'FontSize', base_font_size); 
ylabel('Power (kW)', 'FontSize', base_font_size);
lgd = legend('Day-Ahead Plan', 'Intra-day Actual', 'Location', 'northeast'); % MODIFIED: Location
set(lgd, 'FontSize', base_font_size);
set(gca, 'FontSize', base_font_size);
xlim([0, 24]);
% MODIFIED: Add code to expand Y-axis (Top only)
ax = gca; ylim_vals = ax.YLim;
if ylim_vals(2) < 0.1, ylim_vals(2) = 1; end;
if ylim_vals(1) > -0.1 && ylim_vals(1) < 0.1, ax.YLim = [0, ylim_vals(2) * 1.15]; else, ax.YLim = [ylim_vals(1), ylim_vals(2) * 1.15]; end

figure('Name', 'MG3 Discharging Service to SESO');
stairs(time_axis_da, plan_DA.P_mg3_lease_d, 'r--', 'LineWidth', 1.5); hold on;
plot(time_axis_rt, actual_operation.P_mg3_seso_d, 'b-', 'LineWidth', 1.5);
title('MG3 Discharging Service to SESO', 'FontSize', title_font_size); 
xlabel('Time (h)', 'FontSize', base_font_size); 
ylabel('Power (kW)', 'FontSize', base_font_size);
lgd = legend('Day-Ahead Plan', 'Intra-day Actual', 'Location', 'northeast'); % MODIFIED: Location
set(lgd, 'FontSize', base_font_size);
set(gca, 'FontSize', base_font_size);
xlim([0, 24]);
% MODIFIED: Add code to expand Y-axis (Top only)
ax = gca; ylim_vals = ax.YLim;
if ylim_vals(2) < 0.1, ylim_vals(2) = 1; end;
if ylim_vals(1) > -0.1 && ylim_vals(1) < 0.1, ax.YLim = [0, ylim_vals(2) * 1.15]; else, ax.YLim = [ylim_vals(1), ylim_vals(2) * 1.15]; end

% --- SESO 自有储能状态对比 ---
% --- SESO Own Storage Status ---
figure('Name', 'SESO Own Storage SOC Comparison');
stairs(time_axis_da, plan_DA.E_seso(1:24), 'r--', 'LineWidth', 1.5); hold on;
plot(time_axis_rt, actual_operation.E_seso, 'b-', 'LineWidth', 1.5);
title('SESO Own Storage SOC (DA vs ID)', 'FontSize', title_font_size);
xlabel('Time (h)', 'FontSize', base_font_size); 
ylabel('Energy (kWh)', 'FontSize', base_font_size);
lgd = legend('Day-Ahead Plan', 'Intra-day Actual', 'Location', 'northeast'); % MODIFIED: Location
set(lgd, 'FontSize', base_font_size);
set(gca, 'FontSize', base_font_size);
xlim([0, 24]);
% MODIFIED: Add code to expand Y-axis (Top only)
ax = gca; ylim_vals = ax.YLim;
if ylim_vals(2) < 0.1, ylim_vals(2) = 1; end;
if ylim_vals(1) > -0.1 && ylim_vals(1) < 0.1, ax.YLim = [0, ylim_vals(2) * 1.15]; else, ax.YLim = [ylim_vals(1), ylim_vals(2) * 1.15]; end

figure('Name', 'SESO Own Storage Power Comparison');
stairs(time_axis_da, plan_DA.P_seso_ch, 'r--', 'LineWidth', 1.5); hold on;
plot(time_axis_rt, actual_operation.P_seso_ch, 'b-', 'LineWidth', 1.5);
stairs(time_axis_da, -plan_DA.P_seso_dis, 'm--', 'LineWidth', 1.5);
plot(time_axis_rt, -actual_operation.P_seso_dis, 'c-', 'LineWidth', 1.5);
yline(0, 'k:', 'LineWidth', 1);
title('SESO Own Storage Power (DA vs ID)', 'FontSize', title_font_size);
xlabel('Time (h)', 'FontSize', base_font_size); 
ylabel('Power (kW)', 'FontSize', base_font_size);
lgd = legend('DA Charge', 'ID Charge', 'DA Discharge', 'ID Discharge', 'Location', 'northeast'); % MODIFIED: Location
set(lgd, 'FontSize', base_font_size);
set(gca, 'FontSize', base_font_size);
xlim([0, 24]);
% MODIFIED: Add code to expand Y-axis (Both)
ax = gca; ylim_vals = ax.YLim;
if ylim_vals(1) > -0.1, ylim_vals(1) = -1; end;
if ylim_vals(2) < 0.1, ylim_vals(2) = 1; end;
ax.YLim = [ylim_vals(1) * 1.15, ylim_vals(2) * 1.15];

% --- MG1 自有储能状态对比 ---
% --- MG1 Own Storage Status ---
figure('Name', 'MG1 Own Storage SOC Comparison');
stairs(time_axis_da, plan_DA.E_mg1, 'r--', 'LineWidth', 1.5); hold on;
plot(time_axis_rt, actual_operation.E_mg1, 'b-', 'LineWidth', 1.5);
title('MG1 Own Storage SOC (DA vs ID)', 'FontSize', title_font_size);
xlabel('Time (h)', 'FontSize', base_font_size); 
ylabel('Energy (kWh)', 'FontSize', base_font_size);
lgd = legend('Day-Ahead Plan', 'Intra-day Actual', 'Location', 'northeast'); % MODIFIED: Location
set(lgd, 'FontSize', base_font_size);
set(gca, 'FontSize', base_font_size);
xlim([0, 24]);
% MODIFIED: Add code to expand Y-axis (Top only)
ax = gca; ylim_vals = ax.YLim;
if ylim_vals(2) < 0.1, ylim_vals(2) = 1; end;
if ylim_vals(1) > -0.1 && ylim_vals(1) < 0.1, ax.YLim = [0, ylim_vals(2) * 1.15]; else, ax.YLim = [ylim_vals(1), ylim_vals(2) * 1.15]; end

figure('Name', 'MG1 Own Storage Power Comparison');
stairs(time_axis_da, plan_DA.P_mg1_batc, 'r--', 'LineWidth', 1.5); hold on;
plot(time_axis_rt, actual_operation.P_mg1_batc, 'b-', 'LineWidth', 1.5);
stairs(time_axis_da, -plan_DA.P_mg1_batd, 'm--', 'LineWidth', 1.5);
plot(time_axis_rt, -actual_operation.P_mg1_batd, 'c-', 'LineWidth', 1.5);
yline(0, 'k:', 'LineWidth', 1);
title('MG1 Own Storage Power (DA vs ID)', 'FontSize', title_font_size);
xlabel('Time (h)', 'FontSize', base_font_size); 
ylabel('Power (kW)', 'FontSize', base_font_size);
lgd = legend('DA Charge', 'ID Charge', 'DA Discharge', 'ID Discharge', 'Location', 'northeast'); % MODIFIED: Location
set(lgd, 'FontSize', base_font_size);
set(gca, 'FontSize', base_font_size);
xlim([0, 24]);
% MODIFIED: Add code to expand Y-axis (Both)
ax = gca; ylim_vals = ax.YLim;
if ylim_vals(1) > -0.1, ylim_vals(1) = -1; end;
if ylim_vals(2) < 0.1, ylim_vals(2) = 1; end;
ax.YLim = [ylim_vals(1) * 1.15, ylim_vals(2) * 1.15];

% --- MG2 自有储能状态对比 ---
% --- MG2 Own Storage Status ---
figure('Name', 'MG2 Own Storage SOC Comparison');
stairs(time_axis_da, plan_DA.E_mg2, 'r--', 'LineWidth', 1.5); hold on;
plot(time_axis_rt, actual_operation.E_mg2, 'b-', 'LineWidth', 1.5);
title('MG2 Own Storage SOC (DA vs ID)', 'FontSize', title_font_size);
xlabel('Time (h)', 'FontSize', base_font_size); 
ylabel('Energy (kWh)', 'FontSize', base_font_size);
lgd = legend('Day-Ahead Plan', 'Intra-day Actual', 'Location', 'northeast'); % MODIFIED: Location
set(lgd, 'FontSize', base_font_size);
set(gca, 'FontSize', base_font_size);
xlim([0, 24]);
% MODIFIED: Add code to expand Y-axis (Top only)
ax = gca; ylim_vals = ax.YLim;
if ylim_vals(2) < 0.1, ylim_vals(2) = 1; end;
if ylim_vals(1) > -0.1 && ylim_vals(1) < 0.1, ax.YLim = [0, ylim_vals(2) * 1.15]; else, ax.YLim = [ylim_vals(1), ylim_vals(2) * 1.15]; end

figure('Name', 'MG2 Own Storage Power Comparison');
stairs(time_axis_da, plan_DA.P_mg2_batc, 'r--', 'LineWidth', 1.5); hold on;
plot(time_axis_rt, actual_operation.P_mg2_batc, 'b-', 'LineWidth', 1.5);
stairs(time_axis_da, -plan_DA.P_mg2_batd, 'm--', 'LineWidth', 1.5);
plot(time_axis_rt, -actual_operation.P_mg2_batd, 'c-', 'LineWidth', 1.5);
yline(0, 'k:', 'LineWidth', 1);
title('MG2 Own Storage Power (DA vs ID)', 'FontSize', title_font_size);
xlabel('Time (h)','FontSize', base_font_size); 
ylabel('Power (kW)', 'FontSize', base_font_size);
lgd = legend('DA Charge', 'ID Charge', 'DA Discharge', 'ID Discharge', 'Location', 'northeast'); % MODIFIED: Location
set(lgd, 'FontSize', base_font_size);
set(gca, 'FontSize', base_font_size);
xlim([0, 24]);
% MODIFIED: Add code to expand Y-axis (Both)
ax = gca; ylim_vals = ax.YLim;
if ylim_vals(1) > -0.1, ylim_vals(1) = -1; end;
if ylim_vals(2) < 0.1, ylim_vals(2) = 1; end;
ax.YLim = [ylim_vals(1) * 1.15, ylim_vals(2) * 1.15];

% --- MG3 自有储能状态对比 ---
% --- MG3 Own Storage Status ---
figure('Name', 'MG3 Own Storage SOC Comparison');
stairs(time_axis_da, plan_DA.E_mg3, 'r--', 'LineWidth', 1.5); hold on;
plot(time_axis_rt, actual_operation.E_mg3, 'b-', 'LineWidth', 1.5);
title('MG3 Own Storage SOC (DA vs ID)', 'FontSize', title_font_size);
xlabel('Time (h)', 'FontSize', base_font_size); 
ylabel('Energy (kWh)', 'FontSize', base_font_size);
lgd = legend('Day-Ahead Plan', 'Intra-day Actual', 'Location', 'northeast'); % MODIFIED: Location
set(lgd, 'FontSize', base_font_size);
set(gca, 'FontSize', base_font_size);
xlim([0, 24]);
% MODIFIED: Add code to expand Y-axis (Top only)
ax = gca; ylim_vals = ax.YLim;
if ylim_vals(2) < 0.1, ylim_vals(2) = 1; end;
if ylim_vals(1) > -0.1 && ylim_vals(1) < 0.1, ax.YLim = [0, ylim_vals(2) * 1.15]; else, ax.YLim = [ylim_vals(1), ylim_vals(2) * 1.15]; end

figure('Name', 'MG3 Own Storage Power Comparison');
stairs(time_axis_da, plan_DA.P_mg3_batc, 'r--', 'LineWidth', 1.5); hold on;
plot(time_axis_rt, actual_operation.P_mg3_batc, 'b-', 'LineWidth', 1.5);
stairs(time_axis_da, -plan_DA.P_mg3_batd, 'm--', 'LineWidth', 1.5);
plot(time_axis_rt, -actual_operation.P_mg3_batd, 'c-', 'LineWidth', 1.5);
yline(0, 'k:', 'LineWidth', 1);
title('MG3 Own Storage Power (DA vs ID)', 'FontSize', title_font_size);
xlabel('Time (h)', 'FontSize', base_font_size); 
ylabel('Power (kW)', 'FontSize', base_font_size);
lgd = legend('DA Charge', 'ID Charge', 'DA Discharge', 'ID Discharge', 'Location', 'northeast'); % MODIFIED: Location
set(lgd, 'FontSize', base_font_size);
set(gca, 'FontSize', base_font_size);
xlim([0, 24]);
% MODIFIED: Add code to expand Y-axis (Both)
ax = gca; ylim_vals = ax.YLim;
if ylim_vals(1) > -0.1, ylim_vals(1) = -1; end;
if ylim_vals(2) < 0.1, ylim_vals(2) = 1; end;
ax.YLim = [ylim_vals(1) * 1.15, ylim_vals(2) * 1.15];


fprintf('Program finished.\n');


%% 辅助函数 (保持不变)
function sub_data = filter_data_by_window(full_data, window_idx)
    fields = fieldnames(full_data);
    sub_data = struct();
    for i = 1:length(fields)
        field = fields{i};
        if ismatrix(full_data.(field)) && size(full_data.(field), 2) == 96
             sub_data.(field) = full_data.(field)(:, window_idx);
        elseif isvector(full_data.(field)) && length(full_data.(field)) == 96
             sub_data.(field) = full_data.(field)(window_idx);
        else
            % 如果不是时间序列数据 (例如某些固定参数)，直接复制
             sub_data.(field) = full_data.(field);
        end
    end
end


function zero_deltas = get_zero_deltas(window_len)
    % 生成一个所有调整量都为0的结构体
    zero_deltas.dso_grid   = zeros(1, window_len);
    zero_deltas.dso_seso_c = zeros(1, window_len);
    zero_deltas.dso_seso_d = zeros(1, window_len);
    zero_deltas.mg1_net    = zeros(1, window_len);
    zero_deltas.mg2_net    = zeros(1, window_len);
    zero_deltas.mg3_net    = zeros(1, window_len);
    zero_deltas.mg1_seso_c = zeros(1, window_len);
    zero_deltas.mg1_seso_d = zeros(1, window_len);
    zero_deltas.mg2_seso_c = zeros(1, window_len);
    zero_deltas.mg2_seso_d = zeros(1, window_len);
    zero_deltas.mg3_seso_c = zeros(1, window_len);
    zero_deltas.mg3_seso_d = zeros(1, window_len);
    
    % 在优化失败时，返回的储能状态设为 NaN，由主程序处理
    zero_deltas.E_seso_actual_t1 = NaN;
    zero_deltas.P_seso_ch_actual_t1 = 0;
    zero_deltas.P_seso_dis_actual_t1 = 0;
    zero_deltas.E_mg1_actual_t1 = NaN;
    zero_deltas.P_mg1_batc_actual_t1 = 0;
    zero_deltas.P_mg1_batd_actual_t1 = 0;
    zero_deltas.E_mg2_actual_t1 = NaN;
    zero_deltas.P_mg2_batc_actual_t1 = 0;
    zero_deltas.P_mg2_batd_actual_t1 = 0;
    zero_deltas.E_mg3_actual_t1 = NaN;
    zero_deltas.P_mg3_batc_actual_t1 = 0;
    zero_deltas.P_mg3_batd_actual_t1 = 0;
end
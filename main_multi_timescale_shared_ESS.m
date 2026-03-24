%% =================================================================
%  考虑共享储能的配电网-多微网多时间尺度调度主程序
% =================================================================
clc; clear; close all;
addpath(genpath(pwd));

%% 参数设置
global T H num_steps_RT maxIter_intraday rho_intraday tolerant_intraday;
T = 24;
H = 4; % 滚动优化时域 (1小时)
num_steps_RT = 96; % 15分钟间隔
maxIter_intraday = 20;
rho_intraday = 1;
tolerant_intraday = 1;

%% 1. 准备数据
disp('--- 准备仿真数据 ---');
[forecasts_DA, forecasts_RT] = prepare_data_shared_ESS();

%% 2. 执行日前优化 (包含共享储能租赁)
disp('--- 开始执行日前优化 ---');
[plan_DA, shared_ess_plan, admm_results_DA] = run_day_ahead_shared_ESS(forecasts_DA);
if isempty(plan_DA)
    error('日前优化失败，程序终止。');
end

%% 3. 初始化实际运行结果存储
disp('--- 初始化运行变量 ---');
actual = initialize_actual_variables(num_steps_RT, plan_DA);
violation_costs = zeros(1, num_steps_RT);
convergence_info = [];

%% 4. 日内滚动优化 (分布式ADMM求解)
disp('--- 开始执行日内滚动优化 ---');
fprintf('总时段数: %d, 滚动时域: %d\n', num_steps_RT, H);

for k = 1:num_steps_RT
    fprintf('\n>>> 时段 %d/%d - 正在执行滚动优化...\n', k, num_steps_RT);
    
    % 4.1 提取当前窗口的超短期预测数据
    window_forecasts = get_window_forecasts(forecasts_RT, k, H);
    
    % 4.2 调用分布式日内优化
    [deltas, violation_cost, conv_info] = run_intraday_distributed_ADMM(...
        k, actual, plan_DA, shared_ess_plan, window_forecasts, H, ...
        maxIter_intraday, rho_intraday, tolerant_intraday);
    
    % 4.3 执行并更新系统状态
    if ~isempty(deltas)
        actual = update_actual_state(actual, deltas, plan_DA, shared_ess_plan, k);
        violation_costs(k) = violation_cost;
        convergence_info = [convergence_info; conv_info];
        
        fprintf('    调整完成 - 违约成本: %.2f, 迭代次数: %d\n', ...
                violation_cost, conv_info.iterations);
    else
        warning('时段 %d 优化失败，使用日前计划', k);
        violation_costs(k) = 0;
    end
    
    % 4.4 显示进度
    if mod(k, 16) == 0
        fprintf('已完成 %.1f%%\n', k/num_steps_RT*100);
    end
end

disp('--- 日内滚动优化全部完成 ---');

%% 5. 计算总体性能指标
performance = calculate_performance_metrics(actual, plan_DA, violation_costs, forecasts_RT);

%% 6. 结果可视化
plot_multi_timescale_results(actual, plan_DA, shared_ess_plan, ...
                            violation_costs, convergence_info, performance);

%% 7. 保存结果
save_results(actual, plan_DA, shared_ess_plan, violation_costs, performance);

fprintf('\n======= 调度完成 =======\n');
fprintf('总运行成本: %.2f 元\n', performance.total_cost);
fprintf('违约成本: %.2f 元\n', performance.violation_cost);
fprintf('调整成本: %.2f 元\n', performance.adjustment_cost);
fprintf('平均迭代次数: %.1f\n', performance.avg_iterations);
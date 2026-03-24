%% 主函数：计及共享储能的配电网-多微网多时间尺度协同调度 (一键全流程)
clc; clear; close all;

fprintf('======================================================\n');
fprintf('      配电网-微网多时间尺度协同调度程序 V2.0 (集中式日内) \n');
fprintf('======================================================\n\n');

%% 1. 系统参数定义
T_DA = 1:24;
T_RT = (15/60):(15/60):24; % 96个15分钟点
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
save('plan_DA.mat', 'plan_DA', 'admm_results_DA', 'T_DA', 'T_RT');
fprintf('日前调度完成，已生成24小时运行计划，并保存至 plan_DA.mat。\n\n');

%% 4. 日内滚动调度 (Intra-day Rolling Scheduling)
fprintf('(3/4) 开始执行日内滚动调度 (总共 %d 个时间点)...\n', length(T_RT));
actual_operation = run_intraday_process(plan_DA, dso_data, mg_data, T_RT, H);
fprintf('日内滚动调度全部完成。\n\n');

save('intraday_results.mat', 'actual_operation', 'plan_DA');
fprintf('日内结果已保存到 intraday_results.mat。\n');

fprintf('(4/4) 全流程完成。\n');

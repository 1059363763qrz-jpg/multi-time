%% 单独运行日内滚动调度（直接读取日前结果）
clc; clear; close all;

fprintf('========================================\n');
fprintf('      日内调度独立运行入口\n');
fprintf('========================================\n\n');

if ~isfile('plan_DA.mat')
    error('未找到 plan_DA.mat，请先运行 run_day_ahead_only.m');
end

S = load('plan_DA.mat');
if ~isfield(S, 'plan_DA')
    error('plan_DA.mat 中缺少变量 plan_DA。');
end

plan_DA = S.plan_DA;
if isfield(S, 'admm_results_DA')
    admm_results_DA = S.admm_results_DA; %#ok<NASGU>
else
    admm_results_DA = []; %#ok<NASGU>
end
if isfield(S, 'T_DA')
    T_DA = S.T_DA;
else
    T_DA = 1:24;
end
if isfield(S, 'T_RT')
    T_RT = S.T_RT;
else
    T_RT = (15/60):(15/60):24;
end

H = 4;

fprintf('(1/3) 已加载 plan_DA.mat，正在准备日内数据...\n');
dso_data = prepare_data_DSO(T_DA, T_RT);
mg_data = prepare_data_MGs(T_DA, T_RT);
fprintf('数据准备完成。\n\n');

fprintf('(2/3) 开始执行日内滚动调度 (总共 %d 个时间点)...\n', length(T_RT));
actual_operation = run_intraday_process(plan_DA, dso_data, mg_data, T_RT, H);
fprintf('日内滚动调度全部完成。\n\n');

fprintf('(3/3) 正在保存日内结果...\n');
save('intraday_results.mat', 'actual_operation', 'plan_DA');
fprintf('日内结果已保存到 intraday_results.mat。\n');

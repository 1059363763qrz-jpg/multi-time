%% 单独运行日前调度并保存结果
clc; clear; close all;

fprintf('========================================\n');
fprintf('      日前调度独立运行入口\n');
fprintf('========================================\n\n');

T_DA = 1:24;
T_RT = (15/60):(15/60):24;

fprintf('(1/3) 正在准备数据...\n');
dso_data = prepare_data_DSO(T_DA, T_RT);
mg_data = prepare_data_MGs(T_DA, T_RT);
fprintf('数据准备完成。\n\n');

fprintf('(2/3) 正在执行日前调度...\n');
[plan_DA, admm_results_DA] = run_day_ahead_scheduling(dso_data.forecasts_DA, mg_data.forecasts_DA);
if isempty(plan_DA)
    error('日前调度失败：plan_DA 为空，程序终止。');
end
fprintf('日前调度求解成功。\n\n');

fprintf('(3/3) 正在保存日前结果...\n');
save('plan_DA.mat', 'plan_DA', 'admm_results_DA', 'T_DA', 'T_RT');
fprintf('日前结果已保存到 plan_DA.mat。\n');

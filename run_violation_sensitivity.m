%% 违约惩罚因子灵敏度分析
clc; clear; close all;

fprintf('========================================\n');
fprintf('      违约惩罚因子灵敏度分析\n');
fprintf('========================================\n\n');

check_solver_available();

plan_file = 'plan_DA.mat';
rt_file = 'rt_scenario_fixed.mat';
result_root = fullfile('results', 'violation_sensitivity');
plot_dir = fullfile(result_root, 'plots');

if ~isfile(plan_file)
    error('未找到 plan_DA.mat，请先运行 run_day_ahead_only.m');
end

S = load(plan_file);
if ~isfield(S, 'plan_DA')
    error('plan_DA.mat 中缺少 plan_DA。');
end
plan_DA = S.plan_DA;
if isfield(S, 'T_DA'), T_DA = S.T_DA; else, T_DA = 1:24; end
if isfield(S, 'T_RT'), T_RT = S.T_RT; else, T_RT = (15/60):(15/60):24; end

if ~isfile(rt_file)
    fprintf('未找到固定RT场景，正在创建 %s ...\n', rt_file);
    generate_fixed_rt_scenario(rt_file, 2026, T_DA, T_RT);
end

R = load(rt_file);
required_fields = {'dso_data', 'mg_data'};
for i = 1:numel(required_fields)
    if ~isfield(R, required_fields{i})
        error('固定场景文件 %s 缺少字段 %s。', rt_file, required_fields{i});
    end
end
dso_data = R.dso_data;
mg_data = R.mg_data;
if isfield(R, 'rng_seed'), rng_seed = R.rng_seed; else, rng_seed = NaN; end

violation_list = [0, 1, 5, 10, 20, 50, 100, 200, 500, 1000];
H = 4;

if ~exist(result_root, 'dir'), mkdir(result_root); end
if ~exist(plot_dir, 'dir'), mkdir(plot_dir); end

summary_rows = [];
summary_table = table();

global_tic = tic;
for i = 1:numel(violation_list)
    penalty_value = violation_list(i);
    fprintf('\n==== Penalty = %.4g (%d/%d) ====\n', penalty_value, i, numel(violation_list));

    params = struct('violation', penalty_value, 'fail_fast', true);
    [actual_operation, run_stats] = run_intraday_process(plan_DA, dso_data, mg_data, T_RT, H, params);

    violation_series = run_stats.violation_series;
    summary = build_summary(plan_DA, actual_operation, run_stats, penalty_value);
    summary_rows = [summary_rows; summary]; %#ok<AGROW>

    penalty_dir = fullfile(result_root, sprintf('penalty_%g', penalty_value));
    if ~exist(penalty_dir, 'dir'), mkdir(penalty_dir); end

    save(fullfile(penalty_dir, 'result.mat'), 'penalty_value', 'plan_DA', 'actual_operation', 'violation_series', 'summary', 'run_stats');

    ts_tbl = build_timeseries_table(T_RT, plan_DA, actual_operation, violation_series);
    writetable(ts_tbl, fullfile(penalty_dir, 'timeseries.csv'));

    s_tbl = struct2table(summary);
    writetable(s_tbl, fullfile(penalty_dir, 'summary.csv'));

    summary_table = [summary_table; s_tbl]; %#ok<AGROW>
end

writetable(summary_table, fullfile(result_root, 'summary_all.csv'));
summary_all = summary_table; %#ok<NASGU>
save(fullfile(result_root, 'summary_all.mat'), 'summary_all', 'summary_rows');

metadata_file = fullfile(result_root, 'metadata.txt');
fid = fopen(metadata_file, 'w');
fprintf(fid, 'generated_at=%s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, 'plan_file=%s\n', plan_file);
fprintf(fid, 'rt_file=%s\n', rt_file);
fprintf(fid, 'rng_seed=%s\n', num2str(rng_seed));
fprintf(fid, 'T_DA_len=%d\n', numel(T_DA));
fprintf(fid, 'T_RT_len=%d\n', numel(T_RT));
fprintf(fid, 'violation_list=%s\n', mat2str(violation_list));
fprintf(fid, 'note_total_violation=raw_time_step_sum\n');
fprintf(fid, 'note_total_violation_energy=0.25*total_violation (15-min energy)\n');
fprintf(fid, 'elapsed_sec=%.6f\n', toc(global_tic));
fclose(fid);

save(fullfile(result_root, 'config.mat'), 'plan_file', 'rt_file', 'violation_list', 'T_DA', 'T_RT', 'rng_seed');

try
    make_basic_plots(summary_table, plot_dir);
catch ME
    warning('绘图失败（不影响主流程）：%s', ME.message);
end

fprintf('\n灵敏度分析完成。结果目录：%s\n', result_root);

%% ---------------- 本脚本内部函数 ----------------
function summary = build_summary(plan_DA, actual_operation, run_stats, penalty_value)
    v = run_stats.violation_series;

    mg1_c = nansum(v.v_mg1_c); mg1_d = nansum(v.v_mg1_d);
    mg2_c = nansum(v.v_mg2_c); mg2_d = nansum(v.v_mg2_d);
    mg3_c = nansum(v.v_mg3_c); mg3_d = nansum(v.v_mg3_d);
    seso_c = nansum(v.v_seso_c); seso_d = nansum(v.v_seso_d);

    violation_charge_total = mg1_c + mg2_c + mg3_c + seso_c;
    violation_discharge_total = mg1_d + mg2_d + mg3_d + seso_d;
    total_violation = violation_charge_total + violation_discharge_total;

    % 能量口径(15分钟=0.25h)
    violation_charge_energy = 0.25 * violation_charge_total;
    violation_discharge_energy = 0.25 * violation_discharge_total;
    total_violation_energy = 0.25 * total_violation;

    % 日前承诺总量（功率时序求和口径）
    mg1_commitment_total = nansum(plan_DA.P_mg1_lease_c + plan_DA.P_mg1_lease_d);
    mg2_commitment_total = nansum(plan_DA.P_mg2_lease_c + plan_DA.P_mg2_lease_d);
    mg3_commitment_total = nansum(plan_DA.P_mg3_lease_c + plan_DA.P_mg3_lease_d);
    seso_commitment_total = nansum(plan_DA.P_seso_ch + plan_DA.P_seso_dis);
    total_commitment_mg_lease = mg1_commitment_total + mg2_commitment_total + mg3_commitment_total;
    total_commitment_all = total_commitment_mg_lease + seso_commitment_total;

    mg1_violation_total = mg1_c + mg1_d;
    mg2_violation_total = mg2_c + mg2_d;
    mg3_violation_total = mg3_c + mg3_d;
    seso_violation_total = seso_c + seso_d;

    mg1_violation_energy = 0.25 * mg1_violation_total;
    mg2_violation_energy = 0.25 * mg2_violation_total;
    mg3_violation_energy = 0.25 * mg3_violation_total;
    seso_violation_energy = 0.25 * seso_violation_total;

    total_violation_mg_only = mg1_violation_total + mg2_violation_total + mg3_violation_total;
    total_violation_energy_mg_only = 0.25 * total_violation_mg_only;

    % 安全除法，分母为0时返回0
    violation_rate_all = safe_div(total_violation_energy, total_commitment_all);
    violation_rate_mg_only = safe_div(total_violation_energy_mg_only, total_commitment_mg_lease);
    mg1_violation_rate = safe_div(mg1_violation_energy, mg1_commitment_total);
    mg2_violation_rate = safe_div(mg2_violation_energy, mg2_commitment_total);
    mg3_violation_rate = safe_div(mg3_violation_energy, mg3_commitment_total);
    seso_violation_rate = safe_div(seso_violation_energy, seso_commitment_total);

    penalty_cost = nansum(run_stats.cost_step.penalty_cost);
    adjustment_cost = nansum(run_stats.cost_step.adjustment_cost);
    objective_total = nansum(run_stats.cost_step.objective_total);
    actual_operating_cost_excl_penalty = nansum(run_stats.cost_step.actual_operating_cost_excl_penalty);

    % 资源行为指标
    grid_purchase_total = nansum(max(actual_operation.P_dso_grid, 0));
    grid_purchase_peak = max(actual_operation.P_dso_grid);
    mg1_gt_total = nansum(actual_operation.P_mg1_gt);
    mg2_gt_total = nansum(actual_operation.P_mg2_gt);
    mg3_gt_total = nansum(actual_operation.P_mg3_gt);

    mg1_battery_charge_total = nansum(actual_operation.P_mg1_batc);
    mg1_battery_discharge_total = nansum(actual_operation.P_mg1_batd);
    mg2_battery_charge_total = nansum(actual_operation.P_mg2_batc);
    mg2_battery_discharge_total = nansum(actual_operation.P_mg2_batd);
    mg3_battery_charge_total = nansum(actual_operation.P_mg3_batc);
    mg3_battery_discharge_total = nansum(actual_operation.P_mg3_batd);
    seso_charge_total = nansum(actual_operation.P_seso_ch);
    seso_discharge_total = nansum(actual_operation.P_seso_dis);

    mg1_battery_throughput = mg1_battery_charge_total + mg1_battery_discharge_total;
    mg2_battery_throughput = mg2_battery_charge_total + mg2_battery_discharge_total;
    mg3_battery_throughput = mg3_battery_charge_total + mg3_battery_discharge_total;
    seso_throughput = seso_charge_total + seso_discharge_total;

    % 需求响应指标（来自日内一阶段实际量）
    load_cut_total = nansum(abs(actual_operation.P_mg1_e_cut)) + nansum(abs(actual_operation.P_mg2_e_cut)) + nansum(abs(actual_operation.P_mg3_e_cut));
    load_shift_total = nansum(abs(actual_operation.P_mg1_e_tran)) + nansum(abs(actual_operation.P_mg2_e_tran)) + nansum(abs(actual_operation.P_mg3_e_tran));

    summary = struct();
    summary.penalty_value = penalty_value;
    summary.success = run_stats.success;
    summary.failed_steps_count = run_stats.failed_steps_count;
    summary.solve_time_total = run_stats.solve_time_total;
    summary.solve_time_avg = run_stats.solve_time_avg;

    summary.total_violation = total_violation;
    summary.total_violation_energy = total_violation_energy;
    summary.violation_charge_total = violation_charge_total;
    summary.violation_discharge_total = violation_discharge_total;
    summary.violation_charge_energy = violation_charge_energy;
    summary.violation_discharge_energy = violation_discharge_energy;

    summary.mg1_violation_total = mg1_violation_total;
    summary.mg2_violation_total = mg2_violation_total;
    summary.mg3_violation_total = mg3_violation_total;
    summary.seso_violation_total = seso_violation_total;
    summary.mg1_violation_energy = mg1_violation_energy;
    summary.mg2_violation_energy = mg2_violation_energy;
    summary.mg3_violation_energy = mg3_violation_energy;
    summary.seso_violation_energy = seso_violation_energy;

    summary.total_commitment_all = total_commitment_all;
    summary.total_commitment_mg_lease = total_commitment_mg_lease;
    summary.total_violation_mg_only = total_violation_mg_only;
    summary.total_violation_energy_mg_only = total_violation_energy_mg_only;
    summary.violation_rate_all = violation_rate_all;
    summary.violation_rate_mg_only = violation_rate_mg_only;

    summary.mg1_commitment_total = mg1_commitment_total;
    summary.mg2_commitment_total = mg2_commitment_total;
    summary.mg3_commitment_total = mg3_commitment_total;
    summary.seso_commitment_total = seso_commitment_total;
    summary.mg1_violation_rate = mg1_violation_rate;
    summary.mg2_violation_rate = mg2_violation_rate;
    summary.mg3_violation_rate = mg3_violation_rate;
    summary.seso_violation_rate = seso_violation_rate;

    summary.mg1_charge_violation_total = mg1_c;
    summary.mg1_discharge_violation_total = mg1_d;
    summary.mg2_charge_violation_total = mg2_c;
    summary.mg2_discharge_violation_total = mg2_d;
    summary.mg3_charge_violation_total = mg3_c;
    summary.mg3_discharge_violation_total = mg3_d;
    summary.seso_charge_violation_total = seso_c;
    summary.seso_discharge_violation_total = seso_d;

    summary.penalty_cost = penalty_cost;
    summary.adjustment_cost = adjustment_cost;
    summary.objective_total = objective_total;
    summary.actual_operating_cost_excl_penalty = actual_operating_cost_excl_penalty;

    summary.grid_purchase_total = grid_purchase_total;
    summary.grid_purchase_peak = grid_purchase_peak;
    summary.mg1_gt_total = mg1_gt_total;
    summary.mg2_gt_total = mg2_gt_total;
    summary.mg3_gt_total = mg3_gt_total;

    summary.mg1_battery_charge_total = mg1_battery_charge_total;
    summary.mg1_battery_discharge_total = mg1_battery_discharge_total;
    summary.mg2_battery_charge_total = mg2_battery_charge_total;
    summary.mg2_battery_discharge_total = mg2_battery_discharge_total;
    summary.mg3_battery_charge_total = mg3_battery_charge_total;
    summary.mg3_battery_discharge_total = mg3_battery_discharge_total;
    summary.seso_charge_total = seso_charge_total;
    summary.seso_discharge_total = seso_discharge_total;

    summary.mg1_battery_throughput = mg1_battery_throughput;
    summary.mg2_battery_throughput = mg2_battery_throughput;
    summary.mg3_battery_throughput = mg3_battery_throughput;
    summary.seso_throughput = seso_throughput;

    summary.load_cut_total = load_cut_total;
    summary.load_shift_total = load_shift_total;
end

function ts_tbl = build_timeseries_table(T_RT, plan_DA, actual_operation, violation_series)
    n = numel(T_RT);
    hour_idx = floor((0:n-1)/4) + 1;

    ts_tbl = table();
    ts_tbl.t_idx = (1:n)';
    ts_tbl.time_h = T_RT(:);

    ts_tbl.P_dso_grid_actual = actual_operation.P_dso_grid(:);
    ts_tbl.P_dso_seso_c_actual = actual_operation.P_dso_seso_c(:);
    ts_tbl.P_dso_seso_d_actual = actual_operation.P_dso_seso_d(:);

    ts_tbl.P_mg1_net_actual = actual_operation.P_mg1_net(:);
    ts_tbl.P_mg2_net_actual = actual_operation.P_mg2_net(:);
    ts_tbl.P_mg3_net_actual = actual_operation.P_mg3_net(:);

    ts_tbl.P_mg1_gt_actual = actual_operation.P_mg1_gt(:);
    ts_tbl.P_mg2_gt_actual = actual_operation.P_mg2_gt(:);
    ts_tbl.P_mg3_gt_actual = actual_operation.P_mg3_gt(:);

    ts_tbl.P_mg1_batc_actual = actual_operation.P_mg1_batc(:);
    ts_tbl.P_mg1_batd_actual = actual_operation.P_mg1_batd(:);
    ts_tbl.P_mg2_batc_actual = actual_operation.P_mg2_batc(:);
    ts_tbl.P_mg2_batd_actual = actual_operation.P_mg2_batd(:);
    ts_tbl.P_mg3_batc_actual = actual_operation.P_mg3_batc(:);
    ts_tbl.P_mg3_batd_actual = actual_operation.P_mg3_batd(:);
    ts_tbl.P_seso_ch_actual = actual_operation.P_seso_ch(:);
    ts_tbl.P_seso_dis_actual = actual_operation.P_seso_dis(:);

    ts_tbl.E_mg1_actual = actual_operation.E_mg1(:);
    ts_tbl.E_mg2_actual = actual_operation.E_mg2(:);
    ts_tbl.E_mg3_actual = actual_operation.E_mg3(:);
    ts_tbl.E_seso_actual = actual_operation.E_seso(:);

    ts_tbl.P_mg1_e_cut_actual = actual_operation.P_mg1_e_cut(:);
    ts_tbl.P_mg2_e_cut_actual = actual_operation.P_mg2_e_cut(:);
    ts_tbl.P_mg3_e_cut_actual = actual_operation.P_mg3_e_cut(:);
    ts_tbl.P_mg1_e_tran_actual = actual_operation.P_mg1_e_tran(:);
    ts_tbl.P_mg2_e_tran_actual = actual_operation.P_mg2_e_tran(:);
    ts_tbl.P_mg3_e_tran_actual = actual_operation.P_mg3_e_tran(:);

    ts_tbl.v_mg1_c = violation_series.v_mg1_c(:);
    ts_tbl.v_mg1_d = violation_series.v_mg1_d(:);
    ts_tbl.v_mg2_c = violation_series.v_mg2_c(:);
    ts_tbl.v_mg2_d = violation_series.v_mg2_d(:);
    ts_tbl.v_mg3_c = violation_series.v_mg3_c(:);
    ts_tbl.v_mg3_d = violation_series.v_mg3_d(:);
    ts_tbl.v_seso_c = violation_series.v_seso_c(:);
    ts_tbl.v_seso_d = violation_series.v_seso_d(:);

    ts_tbl.plan_mg1_lease_c = reshape(plan_DA.P_mg1_lease_c(hour_idx), [], 1);
    ts_tbl.plan_mg1_lease_d = reshape(plan_DA.P_mg1_lease_d(hour_idx), [], 1);
    ts_tbl.plan_mg2_lease_c = reshape(plan_DA.P_mg2_lease_c(hour_idx), [], 1);
    ts_tbl.plan_mg2_lease_d = reshape(plan_DA.P_mg2_lease_d(hour_idx), [], 1);
    ts_tbl.plan_mg3_lease_c = reshape(plan_DA.P_mg3_lease_c(hour_idx), [], 1);
    ts_tbl.plan_mg3_lease_d = reshape(plan_DA.P_mg3_lease_d(hour_idx), [], 1);
    ts_tbl.plan_seso_ch = reshape(plan_DA.P_seso_ch(hour_idx), [], 1);
    ts_tbl.plan_seso_dis = reshape(plan_DA.P_seso_dis(hour_idx), [], 1);
end


function y = safe_div(num, den)
    if isempty(den) || isnan(den) || abs(den) < 1e-12
        y = 0;
    else
        y = num / den;
    end
end

function make_basic_plots(summary_table, plot_dir)
    if isempty(summary_table), return; end
    p = summary_table.penalty_value;

    fig = figure('Visible', 'off');
    semilogx(p, summary_table.total_violation, 'o-', 'LineWidth', 1.5);
    xlabel('Penalty'); ylabel('Total violation'); title('Penalty vs Total Violation'); grid on;
    saveas(fig, fullfile(plot_dir, 'penalty_vs_total_violation.png')); close(fig);

    fig = figure('Visible', 'off');
    semilogx(p, summary_table.actual_operating_cost_excl_penalty, 'o-', 'LineWidth', 1.5);
    xlabel('Penalty'); ylabel('Actual operating cost (excl. penalty)'); title('Penalty vs Actual Operating Cost'); grid on;
    saveas(fig, fullfile(plot_dir, 'penalty_vs_actual_operating_cost.png')); close(fig);

    fig = figure('Visible', 'off');
    semilogx(p, summary_table.penalty_cost, 'o-', 'LineWidth', 1.5);
    xlabel('Penalty'); ylabel('Penalty cost'); title('Penalty vs Penalty Cost'); grid on;
    saveas(fig, fullfile(plot_dir, 'penalty_vs_penalty_cost.png')); close(fig);

    fig = figure('Visible', 'off');
    semilogx(p, summary_table.seso_throughput, 'o-', 'LineWidth', 1.5); hold on;
    semilogx(p, summary_table.mg1_battery_throughput, 's-', 'LineWidth', 1.2);
    semilogx(p, summary_table.mg2_battery_throughput, '^-', 'LineWidth', 1.2);
    semilogx(p, summary_table.mg3_battery_throughput, 'd-', 'LineWidth', 1.2);
    xlabel('Penalty'); ylabel('Throughput'); title('Penalty vs Resource Usage'); legend('SESO', 'MG1', 'MG2', 'MG3', 'Location', 'best'); grid on;
    saveas(fig, fullfile(plot_dir, 'penalty_vs_resource_usage.png')); close(fig);

    fig = figure('Visible', 'off');
    scatter(summary_table.total_violation, summary_table.actual_operating_cost_excl_penalty, 60, summary_table.penalty_value, 'filled');
    xlabel('Total violation'); ylabel('Actual operating cost (excl. penalty)'); title('Pareto: Cost vs Violation'); grid on; colorbar;
    saveas(fig, fullfile(plot_dir, 'pareto_cost_vs_violation.png')); close(fig);
end

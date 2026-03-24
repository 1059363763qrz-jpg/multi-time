%% 主函数：计及共享储能的配电网-多微网协同调度 
clc; clear; close all;

%% ADMM 迭代参数
maxIter = 200;      % 最大迭代次数
rho = 1;         % 罚因子 (重要！可能需要根据收敛情况调整)
tolerant = 2;    % 收敛精度
iter = 1;           % 迭代次数初始化

% %% 初始化 primal 和 dual 变量
T = 24;
% === 自适应rho的参数 ===
% mu = 5;      % 原始残差与对偶残差的平衡阈值
% tau_incr = 2; % rho 增大系数
% tau_decr = 2; % rho 减小系数
% rho_history = rho; % 记录rho的变化

% ===  (交互功率) ===
P_dso_charge = zeros(maxIter+1, T); P_dso_discharge = zeros(maxIter+1, T);
P_seso_to_dso_c = zeros(maxIter+1, T); P_seso_to_dso_d = zeros(maxIter+1, T);
P_mg1_lease_c = zeros(maxIter+1, T); P_mg1_lease_d = zeros(maxIter+1, T);
P_mg2_lease_c = zeros(maxIter+1, T); P_mg2_lease_d = zeros(maxIter+1, T);
P_mg3_lease_c = zeros(maxIter+1, T); P_mg3_lease_d = zeros(maxIter+1, T);
P_seso_from_mg1_c = zeros(maxIter+1, T); P_seso_from_mg1_d = zeros(maxIter+1, T);
P_seso_from_mg2_c = zeros(maxIter+1, T); P_seso_from_mg2_d = zeros(maxIter+1, T);
P_seso_from_mg3_c = zeros(maxIter+1, T); P_seso_from_mg3_d = zeros(maxIter+1, T);
P_mg1_net_exchange = zeros(maxIter+1, T);
P_mg2_net_exchange = zeros(maxIter+1, T);
P_mg3_net_exchange = zeros(maxIter+1, T);

% ===  (拉格朗日乘子) ===
lambda_dso_c = zeros(1, T); lambda_dso_d = zeros(1, T);
lambda_mg1_c = zeros(1, T); lambda_mg1_d = zeros(1, T);
lambda_mg2_c = zeros(1, T); lambda_mg2_d = zeros(1, T);
lambda_mg3_c = zeros(1, T); lambda_mg3_d = zeros(1, T);


% 历史数据记录
Obj_DSO = []; Obj_MG1 = []; Obj_MG2 = []; Obj_MG3 = []; Obj_SESO = [];
Residual = [];

%% ADMM 迭代主循环
while iter <= maxIter
    fprintf('正在进行第 %d 次迭代...  %f\n', iter);
    % 步骤 1: 求解微网子问题 (MG1, MG2, MG3)
    [P_mg1_net_exchange(iter+1,:), P_mg1_lease_c(iter+1,:), P_mg1_lease_d(iter+1,:), obj_mg1] = Fun_MG1_Lease(...
        P_seso_from_mg1_c(iter,:), P_seso_from_mg1_d(iter,:), lambda_mg1_c, lambda_mg1_d, rho);
    Obj_MG1(iter) = obj_mg1;

    [P_mg2_net_exchange(iter+1,:), P_mg2_lease_c(iter+1,:), P_mg2_lease_d(iter+1,:), obj_mg2] = Fun_MG2_Lease(...
        P_seso_from_mg2_c(iter,:), P_seso_from_mg2_d(iter,:), lambda_mg2_c, lambda_mg2_d, rho);
    Obj_MG2(iter) = obj_mg2;

    [P_mg3_net_exchange(iter+1,:), P_mg3_lease_c(iter+1,:), P_mg3_lease_d(iter+1,:), obj_mg3] = Fun_MG3_Lease(...
        P_seso_from_mg3_c(iter,:), P_seso_from_mg3_d(iter,:), lambda_mg3_c, lambda_mg3_d, rho);
    Obj_MG3(iter) = obj_mg3;

    % 步骤 2: 求解DSO子问题
    [P_dso_charge(iter+1,:), P_dso_discharge(iter+1,:), obj_dso] = Fun_DSO(...
        P_seso_to_dso_c(iter,:), P_seso_to_dso_d(iter,:), lambda_dso_c, lambda_dso_d, ...
        P_mg1_net_exchange(iter+1,:), P_mg2_net_exchange(iter+1,:), P_mg3_net_exchange(iter+1,:), rho);
    Obj_DSO(iter) = obj_dso;

    % 步骤 3: 求解SESO子问题
    [P_seso_to_dso_c(iter+1,:), P_seso_to_dso_d(iter+1,:), ...
     P_seso_from_mg1_c(iter+1,:), P_seso_from_mg1_d(iter+1,:), ...
     P_seso_from_mg2_c(iter+1,:), P_seso_from_mg2_d(iter+1,:), ...
     P_seso_from_mg3_c(iter+1,:), P_seso_from_mg3_d(iter+1,:), ...
     obj_seso] = Fun_SESO(...
        P_dso_charge(iter+1,:), P_dso_discharge(iter+1,:), lambda_dso_c, lambda_dso_d, ...
        P_mg1_lease_c(iter+1,:), P_mg1_lease_d(iter+1,:), lambda_mg1_c, lambda_mg1_d, ...
        P_mg2_lease_c(iter+1,:), P_mg2_lease_d(iter+1,:), lambda_mg2_c, lambda_mg2_d, ...
        P_mg3_lease_c(iter+1,:), P_mg3_lease_d(iter+1,:), lambda_mg3_c, lambda_mg3_d, rho);
    Obj_SESO(iter) = obj_seso;

    % 步骤 4: 更新拉格朗日乘子
    lambda_dso_c = lambda_dso_c + rho * (P_dso_charge(iter+1,:) - P_seso_to_dso_c(iter+1,:));
    lambda_dso_d = lambda_dso_d + rho * (P_dso_discharge(iter+1,:) - P_seso_to_dso_d(iter+1,:));
    lambda_mg1_c = lambda_mg1_c + rho * (P_mg1_lease_c(iter+1,:) - P_seso_from_mg1_c(iter+1,:));
    lambda_mg1_d = lambda_mg1_d + rho * (P_mg1_lease_d(iter+1,:) - P_seso_from_mg1_d(iter+1,:));
    lambda_mg2_c = lambda_mg2_c + rho * (P_mg2_lease_c(iter+1,:) - P_seso_from_mg2_c(iter+1,:));
    lambda_mg2_d = lambda_mg2_d + rho * (P_mg2_lease_d(iter+1,:) - P_seso_from_mg2_d(iter+1,:));
    lambda_mg3_c = lambda_mg3_c + rho * (P_mg3_lease_c(iter+1,:) - P_seso_from_mg3_c(iter+1,:));
    lambda_mg3_d = lambda_mg3_d + rho * (P_mg3_lease_d(iter+1,:) - P_seso_from_mg3_d(iter+1,:));

    % 步骤 5: 计算并检查收敛性
    res_dso_seso = norm(P_dso_charge(iter+1,:) - P_seso_to_dso_c(iter+1,:)) + norm(P_dso_discharge(iter+1,:) - P_seso_to_dso_d(iter+1,:));
    res_mg1_seso = norm(P_mg1_lease_c(iter+1,:) - P_seso_from_mg1_c(iter+1,:)) + norm(P_mg1_lease_d(iter+1,:) - P_seso_from_mg1_d(iter+1,:));
    res_mg2_seso = norm(P_mg2_lease_c(iter+1,:) - P_seso_from_mg2_c(iter+1,:)) + norm(P_mg2_lease_d(iter+1,:) - P_seso_from_mg2_d(iter+1,:));
    res_mg3_seso = norm(P_mg3_lease_c(iter+1,:) - P_seso_from_mg3_c(iter+1,:)) + norm(P_mg3_lease_d(iter+1,:) - P_seso_from_mg3_d(iter+1,:));

    current_residual = res_dso_seso + res_mg1_seso + res_mg2_seso + res_mg3_seso;
    Residual(iter) = current_residual;
    fprintf('当前残差: %f\n', current_residual);

    if current_residual < tolerant
        fprintf('迭代在第 %d 次收敛！\n', iter);
        break;
    end

    if iter == maxIter
        disp('已达到最大迭代次数，算法未收敛。');
    end

    iter = iter + 1;
end
%% =====================================================================
%                  新增：结果可视化代码
%  =====================================================================

final_iter = iter-1; % 保存最终迭代次数
time_axis = 1:T;
%考虑自适应调节惩罚函数rho
% while iter <= maxIter
%     fprintf('正在进行第 %d 次迭代... rho = %f\n', iter, rho);
%     %步骤 1: 求解微网子问题 (MG1, MG2, MG3)
%     [P_mg1_net_exchange(iter+1,:), P_mg1_lease_c(iter+1,:), P_mg1_lease_d(iter+1,:), obj_mg1] = Fun_MG1_Lease(...
%         P_seso_from_mg1_c(iter,:), P_seso_from_mg1_d(iter,:), lambda_mg1_c, lambda_mg1_d, rho);
%     Obj_MG1(iter) = obj_mg1;
% 
%     [P_mg2_net_exchange(iter+1,:), P_mg2_lease_c(iter+1,:), P_mg2_lease_d(iter+1,:), obj_mg2] = Fun_MG2_Lease(...
%         P_seso_from_mg2_c(iter,:), P_seso_from_mg2_d(iter,:), lambda_mg2_c, lambda_mg2_d, rho);
%     Obj_MG2(iter) = obj_mg2;
% 
%     [P_mg3_net_exchange(iter+1,:), P_mg3_lease_c(iter+1,:), P_mg3_lease_d(iter+1,:), obj_mg3] = Fun_MG3_Lease(...
%         P_seso_from_mg3_c(iter,:), P_seso_from_mg3_d(iter,:), lambda_mg3_c, lambda_mg3_d, rho);
%     Obj_MG3(iter) = obj_mg3;
% 
%     % 步骤 2: 求解DSO子问题
%     [P_dso_charge(iter+1,:), P_dso_discharge(iter+1,:), obj_dso] = Fun_DSO(...
%         P_seso_to_dso_c(iter,:), P_seso_to_dso_d(iter,:), lambda_dso_c, lambda_dso_d, ...
%         P_mg1_net_exchange(iter+1,:), P_mg2_net_exchange(iter+1,:), P_mg3_net_exchange(iter+1,:), rho);
%     Obj_DSO(iter) = obj_dso;
% 
%     %步骤 3: 求解SESO子问题
%     [P_seso_to_dso_c(iter+1,:), P_seso_to_dso_d(iter+1,:), ...
%      P_seso_from_mg1_c(iter+1,:), P_seso_from_mg1_d(iter+1,:), ...
%      P_seso_from_mg2_c(iter+1,:), P_seso_from_mg2_d(iter+1,:), ...
%      P_seso_from_mg3_c(iter+1,:), P_seso_from_mg3_d(iter+1,:), ...
%      obj_seso] = Fun_SESO(...
%         P_dso_charge(iter+1,:), P_dso_discharge(iter+1,:), lambda_dso_c, lambda_dso_d, ...
%         P_mg1_lease_c(iter+1,:), P_mg1_lease_d(iter+1,:), lambda_mg1_c, lambda_mg1_d, ...
%         P_mg2_lease_c(iter+1,:), P_mg2_lease_d(iter+1,:), lambda_mg2_c, lambda_mg2_d, ...
%         P_mg3_lease_c(iter+1,:), P_mg3_lease_d(iter+1,:), lambda_mg3_c, lambda_mg3_d, rho);
%     Obj_SESO(iter) = obj_seso;
% 
%     %步骤 4: 更新拉格朗日乘子
%     lambda_dso_c = lambda_dso_c + rho * (P_dso_charge(iter+1,:) - P_seso_to_dso_c(iter+1,:));
%     lambda_dso_d = lambda_dso_d + rho * (P_dso_discharge(iter+1,:) - P_seso_to_dso_d(iter+1,:));
%     lambda_mg1_c = lambda_mg1_c + rho * (P_mg1_lease_c(iter+1,:) - P_seso_from_mg1_c(iter+1,:));
%     lambda_mg1_d = lambda_mg1_d + rho * (P_mg1_lease_d(iter+1,:) - P_seso_from_mg1_d(iter+1,:));
%     lambda_mg2_c = lambda_mg2_c + rho * (P_mg2_lease_c(iter+1,:) - P_seso_from_mg2_c(iter+1,:));
%     lambda_mg2_d = lambda_mg2_d + rho * (P_mg2_lease_d(iter+1,:) - P_seso_from_mg2_d(iter+1,:));
%     lambda_mg3_c = lambda_mg3_c + rho * (P_mg3_lease_c(iter+1,:) - P_seso_from_mg3_c(iter+1,:));
%     lambda_mg3_d = lambda_mg3_d + rho * (P_mg3_lease_d(iter+1,:) - P_seso_from_mg3_d(iter+1,:));
% 
%     % 步骤 5: 计算原始残差和对偶残差
%     % 原始残差r (Primal Residual)，衡量的是不同主体决策变量间的一致性，这个越大，应该增大惩罚系数
%     res_dso_c_r = P_dso_charge(iter+1,:) - P_seso_to_dso_c(iter+1,:);
%     res_dso_d_r = P_dso_discharge(iter+1,:) - P_seso_to_dso_d(iter+1,:);
%     res_mg1_c_r = P_mg1_lease_c(iter+1,:) - P_seso_from_mg1_c(iter+1,:);
%     res_mg1_d_r = P_mg1_lease_d(iter+1,:) - P_seso_from_mg1_d(iter+1,:);
%     res_mg2_c_r = P_mg2_lease_c(iter+1,:) - P_seso_from_mg2_c(iter+1,:);
%     res_mg2_d_r = P_mg2_lease_d(iter+1,:) - P_seso_from_mg2_d(iter+1,:);
%     res_mg3_c_r = P_mg3_lease_c(iter+1,:) - P_seso_from_mg3_c(iter+1,:);
%     res_mg3_d_r = P_mg3_lease_d(iter+1,:) - P_seso_from_mg3_d(iter+1,:);
% 
%     primal_residual_norm = norm([res_dso_c_r, res_dso_d_r, res_mg1_c_r, res_mg1_d_r, res_mg2_c_r, res_mg2_d_r, res_mg3_c_r, res_mg3_d_r]);
%     Residual(iter) = primal_residual_norm;
%      % 对偶残差s (Dual Residual)，衡量的是是否达到了最优，这个大应该减小惩罚系数
%     dual_res_dso_c_s = rho* (P_seso_to_dso_c(iter+1,:) - P_seso_to_dso_c(iter,:));
%     dual_res_dso_d_s =  rho*(P_seso_to_dso_d(iter+1,:) - P_seso_to_dso_d(iter,:));
%     dual_res_seso_c_s = rho*(P_dso_charge(iter+1,:)-P_dso_charge(iter,:));
%     dual_res_seso_d_s = rho*(P_dso_discharge(iter+1,:)-P_dso_discharge(iter,:));
%     dual_res_mg1_c_s =  rho*(P_seso_from_mg1_c(iter+1,:) - P_seso_from_mg1_c(iter,:));
%     dual_res_mg1_d_s = rho*(P_seso_from_mg1_d(iter+1,:) - P_seso_from_mg1_d(iter,:));
%     dual_res_mg2_c_s =  rho*(P_seso_from_mg2_c(iter+1,:) - P_seso_from_mg2_c(iter,:));
%     dual_res_mg2_d_s =  rho*(P_seso_from_mg2_d(iter+1,:) - P_seso_from_mg2_d(iter,:));
%     dual_res_mg3_c_s = rho*(P_seso_from_mg3_c(iter+1,:) - P_seso_from_mg3_c(iter,:));
%     dual_res_mg3_d_s = rho*(P_seso_from_mg3_d(iter+1,:) - P_seso_from_mg3_d(iter,:));
%  dual_residual_norm = norm([dual_res_dso_c_s, dual_res_dso_d_s,dual_res_seso_c_s, dual_res_seso_d_s,dual_res_mg1_c_s, dual_res_mg1_d_s, dual_res_mg2_c_s,dual_res_mg2_d_s,dual_res_mg3_c_s,dual_res_mg3_d_s]);    
% fprintf('原始残差: %f, 对偶残差: %f\n', primal_residual_norm, dual_residual_norm);
%  % 步骤 6: 检查收敛性
%     if primal_residual_norm < tolerant && dual_residual_norm < tolerant
%         fprintf('迭代在第 %d 次收敛！\n', iter);
%         break;
%     end
% 
%     % 步骤 7: 更新rho
%     mu=5;
%     if primal_residual_norm > mu * dual_residual_norm
%         rho = rho * tau_incr;
%     elseif dual_residual_norm > mu * primal_residual_norm
%         rho = rho / tau_decr;
%     end
%     rho_history(iter+1) = rho; % 记录rho
% 
%     if iter == maxIter
%         disp('已达到最大迭代次数，算法未收敛。');
%     end
% 
%     iter = iter + 1;
% end
% 1. 绘制各主体目标函数值的收敛曲线
%% =====================================================================
%% =====================================================================
%                  Result Visualization (English Version)
%  =====================================================================

final_iter = iter-1; % Store the final iteration number
time_axis = 1:T;

% 1. Plot convergence of objective functions for each agent
figure;
plot(1:final_iter, Obj_DSO(1:final_iter), 'r-o', 'LineWidth', 1.5);
hold on;
plot(1:final_iter, Obj_MG1(1:final_iter), 'g-s', 'LineWidth', 1.5);
plot(1:final_iter, Obj_MG2(1:final_iter), 'b-^', 'LineWidth', 1.5);
plot(1:final_iter, Obj_MG3(1:final_iter), 'm-d', 'LineWidth', 1.5);
plot(1:final_iter, Obj_SESO(1:final_iter), 'k-*', 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Objective Function Value');
title('Convergence of Objective Functions for Each Agent');
legend('DSO', 'Microgrid 1', 'Microgrid 2', 'Microgrid 3', 'SESO', 'Location', 'best');
grid on;
box off;

% 2. Plot convergence of the residual
figure;
plot(1:final_iter, Residual(1:final_iter), 'b-x', 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Residual');
title('ADMM Algorithm Residual Convergence');
grid on;
box off;

% 3. Plot energy storage interaction between DSO and SESO
figure;
bar(time_axis, [P_dso_charge(final_iter+1,:)', -P_dso_discharge(final_iter+1,:)'], 'stacked');
xlabel('Time (h)');
ylabel('Power (kW)');
title('Energy Storage Interaction Power between DSO and SESO');
legend('DSO Charging ', 'DSO Discharging ', 'Location', 'northwest');
grid on;
xlim([0.5, 24.5]);

% 4. Plot energy storage leasing results between Microgrids and SESO
figure;

% Microgrid 1
subplot(3,1,1);
bar(time_axis, [P_mg1_lease_c(final_iter+1,:)', P_mg1_lease_d(final_iter+1,:)']);
ylabel('Power (kW)');
title('Storage Power Leased from Microgrid 1 to SESO');
legend('Charging Service', 'Discharging Service');
grid on;
xlim([0.5, 24.5]);

% Microgrid 2
subplot(3,1,2);
bar(time_axis, [P_mg2_lease_c(final_iter+1,:)', P_mg2_lease_d(final_iter+1,:)']);
ylabel('Power (kW)');
title('Storage Power Leased from Microgrid 2 to SESO');
legend('Charging Service', 'Discharging Service');
grid on;
xlim([0.5, 24.5]);

% Microgrid 3
subplot(3,1,3);
bar(time_axis, [P_mg3_lease_c(final_iter+1,:)', P_mg3_lease_d(final_iter+1,:)']);
xlabel('Time (h)');
ylabel('Power (kW)');
title('Storage Power Leased from Microgrid 3 to SESO');
legend('Charging Service', 'Discharging Service');
grid on;
xlim([0.5, 24.5]);

% 5. Plot net exchange power between microgrids and the distribution grid
figure;
plot(time_axis, P_mg1_net_exchange(final_iter+1,:), 'r-o', 'LineWidth', 1.5);
hold on;
plot(time_axis, P_mg2_net_exchange(final_iter+1,:), 'g-s', 'LineWidth', 1.5);
plot(time_axis, P_mg3_net_exchange(final_iter+1,:), 'b-^', 'LineWidth', 1.5);
yline(0, 'k--', 'LineWidth', 1); % Add a zero line
xlabel('Time (h)');
ylabel('Net Exchange Power (kW)');
title('Net Exchange Power (Positive: Sell, Negative: Buy)');
legend('Microgrid 1', 'Microgrid 2', 'Microgrid 3', 'Location', 'best');
grid on;
xlim([1, 24]);

%% =====================================================================
%                  Detailed Visualization (English Version)
%  =====================================================================

% --- Step 1: Rerun subproblems with converged variables to get detailed results ---
fprintf('Generating final visualization results...\n');

% Get detailed results for MG1
[~, ~, ~, ~, E_bat_mg1, P_batc_mg1, P_batd_mg1, P_e_GT_mg1, P_e_wd_mg1, P_buy_mg1, P_sell_mg1, L_e_mg1] = Fun_MG1_Lease(...
    P_seso_from_mg1_c(final_iter+1,:), P_seso_from_mg1_d(final_iter+1,:), lambda_mg1_c, lambda_mg1_d, rho);

% Get detailed results for MG2
[~, ~, ~, ~, E_bat_mg2, P_batc_mg2, P_batd_mg2, P_e_GT_mg2, P_e_pv_mg2, P_buy_mg2, P_sell_mg2, L_e_mg2] = Fun_MG2_Lease(...
    P_seso_from_mg2_c(final_iter+1,:), P_seso_from_mg2_d(final_iter+1,:), lambda_mg2_c, lambda_mg2_d, rho);

% Get detailed results for MG3
[~, ~, ~, ~, E_bat_mg3, P_batc_mg3, P_batd_mg3, P_e_GT_mg3, P_e_pv_mg3, P_buy_mg3, P_sell_mg3, L_e_mg3] = Fun_MG3_Lease(...
    P_seso_from_mg3_c(final_iter+1,:), P_seso_from_mg3_d(final_iter+1,:), lambda_mg3_c, lambda_mg3_d, rho);

% Get detailed results for DSO
[~, ~, ~, P_grid_exchange, PLOAD_dso, P_loss_dso, P_pv_dso, P_wt_dso, PLOAD_original, Lshift_effect, SIL_effect] = Fun_DSO(...
    P_seso_to_dso_c(final_iter+1,:), P_seso_to_dso_d(final_iter+1,:), lambda_dso_c, lambda_dso_d, ...
    P_mg1_net_exchange(final_iter+1,:), P_mg2_net_exchange(final_iter+1,:), P_mg3_net_exchange(final_iter+1,:), rho);
P_buy_dso = max(0, P_grid_exchange);
P_sell_dso = -min(0, P_grid_exchange);

% Get detailed results for SESO
[~, ~, ~, ~, ~, ~, ~, ~, ~, E_seso, P_seso_ch, P_seso_dis] = Fun_SESO(...
    P_dso_charge(final_iter+1,:), P_dso_discharge(final_iter+1,:), lambda_dso_c, lambda_dso_d, ...
    P_mg1_lease_c(final_iter+1,:), P_mg1_lease_d(final_iter+1,:), lambda_mg1_c, lambda_mg1_d, ...
    P_mg2_lease_c(final_iter+1,:), P_mg2_lease_d(final_iter+1,:), lambda_mg2_c, lambda_mg2_d, ...
    P_mg3_lease_c(final_iter+1,:), P_mg3_lease_d(final_iter+1,:), lambda_mg3_c, lambda_mg3_d, rho);

% --- Step 2: Plot Energy Storage Status ---
figure;
sgtitle('Energy Storage System Status of Each Agent (Converged Result)');

% SESO Own Storage SOC
subplot(4,2,1);
plot(time_axis, E_seso(1:T), 'r-s', 'LineWidth', 1.5);
ylabel('Energy (kWh)');
title('SESO Own Storage SOC');
grid on; xlim([1, 24]);

% SESO Own Storage Charge/Discharge Power
subplot(4,2,2);
bar(time_axis, [P_seso_ch', -P_seso_dis'], 'stacked');
ylabel('Power (kW)');
title('SESO Own Storage Charge/Discharge Power');
legend('Charge', 'Discharge');
grid on; xlim([0.5, 24.5]);

% Microgrid 1 Storage SOC
subplot(4,2,3);
plot(time_axis, E_bat_mg1, 'g-^', 'LineWidth', 1.5);
ylabel('Energy (kWh)');
title('Microgrid 1 Storage SOC');
grid on; xlim([1, 24]);

% Microgrid 1 Storage Charge/Discharge Power
subplot(4,2,4);
bar(time_axis, [P_batc_mg1', -P_batd_mg1'], 'stacked');
ylabel('Power (kW)');
title('Microgrid 1 Storage Charge/Discharge Power');
legend('Charge', 'Discharge');
grid on; xlim([0.5, 24.5]);

% Microgrid 2 Storage SOC
subplot(4,2,5);
plot(time_axis, E_bat_mg2, 'b-o', 'LineWidth', 1.5);
ylabel('Energy (kWh)');
title('Microgrid 2 Storage SOC');
grid on; xlim([1, 24]);

% Microgrid 2 Storage Charge/Discharge Power
subplot(4,2,6);
bar(time_axis, [P_batc_mg2', -P_batd_mg2'], 'stacked');
ylabel('Power (kW)');
title('Microgrid 2 Storage Charge/Discharge Power');
legend('Charge', 'Discharge');
grid on; xlim([0.5, 24.5]);

% Microgrid 3 Storage SOC
subplot(4,2,7);
plot(time_axis, E_bat_mg3, 'm-d', 'LineWidth', 1.5);
xlabel('Time (h)');
ylabel('Energy (kWh)');
title('Microgrid 3 Storage SOC');
grid on; xlim([1, 24]);

% Microgrid 3 Storage Charge/Discharge Power
subplot(4,2,8);
bar(time_axis, [P_batc_mg3', -P_batd_mg3'], 'stacked');
xlabel('Time (h)');
ylabel('Power (kW)');
title('Microgrid 3 Storage Charge/Discharge Power');
legend('Charge', 'Discharge');
grid on; xlim([0.5, 24.5]);

% --- Step 3: Plot Power Balance ---
figure;
sgtitle('Power Balance of Each Agent (Source > 0, Load < 0)');

% DSO Power Balance
subplot(3,2,1);
mg1_source = max(0, P_mg1_net_exchange(final_iter+1,:));
mg1_load = -min(0, P_mg1_net_exchange(final_iter+1,:));
mg2_source = max(0, P_mg2_net_exchange(final_iter+1,:));
mg2_load = -min(0, P_mg2_net_exchange(final_iter+1,:));
mg3_source = max(0, P_mg3_net_exchange(final_iter+1,:));
mg3_load = -min(0, P_mg3_net_exchange(final_iter+1,:));
P_buy_dso = max(0, P_grid_exchange);
P_sell_dso = -min(0, P_grid_exchange);
sources_dso = [P_buy_dso', sum(P_pv_dso,1)', sum(P_wt_dso,1)', mg1_source', mg2_source', mg3_source', P_dso_discharge(final_iter+1,:)'];
loads_dso = [-sum(PLOAD_original,1)', sum(Lshift_effect,1)', sum(SIL_effect,1)', ...
             -P_sell_dso', -P_loss_dso', -mg1_load', -mg2_load', -mg3_load', -P_dso_charge(final_iter+1,:)'];
bar(time_axis, [sources_dso, loads_dso], 'stacked');
title('DSO Power Balance (with DR Decomposition)');
ylabel('Power (kW)');
legend('Grid Purchase','PV','Wind','Buy from MG','SESO Discharge', ...
       'Original Load','Load Shifting (Decrease)','Load Curtailment', ...
       'Sell to Grid','Grid Loss','Sell to MG','SESO Charge', ...
       'Location','bestoutside');
grid on; xlim([0.5, 24.5]);

% SESO Power Balance
subplot(3,2,2);
sources_seso = [P_seso_dis', P_seso_from_mg1_d(final_iter+1,:)', P_seso_from_mg2_d(final_iter+1,:)', P_seso_from_mg3_d(final_iter+1,:)'];
loads_seso = [-P_seso_ch', -P_seso_from_mg1_c(final_iter+1,:)', -P_seso_from_mg2_c(final_iter+1,:)', -P_seso_from_mg3_c(final_iter+1,:)'];
bar(time_axis, [sources_seso, loads_seso], 'stacked');
title('SESO Power Balance');
ylabel('Power (kW)');
legend('Own Dis.','Lease Dis. MG1','Lease Dis. MG2','Lease Dis. MG3', ...
       'Own Ch.','Lease Ch. MG1','Lease Ch. MG2','Lease Ch. MG3','Location','bestoutside');
grid on; xlim([0.5, 24.5]);


% Microgrid 1 Power Balance
subplot(3,2,3);
sources_mg1 = [P_e_GT_mg1', P_e_wd_mg1', P_buy_mg1', P_batd_mg1'];
loads_mg1 = [-L_e_mg1', -P_sell_mg1', -P_batc_mg1'];
bar(time_axis, [sources_mg1, loads_mg1], 'stacked');
title('Microgrid 1 Power Balance');
ylabel('Power (kW)');
legend('Gas Turbine','Wind Power','Grid Purchase','Storage Dis.','Elec. Load','Grid Sale','Storage Ch.','Location','bestoutside');
grid on; xlim([0.5, 24.5]);

% Microgrid 2 Power Balance
subplot(3,2,4);
sources_mg2 = [P_e_GT_mg2', P_e_pv_mg2', P_buy_mg2', P_batd_mg2'];
loads_mg2 = [-L_e_mg2', -P_sell_mg2', -P_batc_mg2'];
bar(time_axis, [sources_mg2, loads_mg2], 'stacked');
title('Microgrid 2 Power Balance');
ylabel('Power (kW)');
legend('Gas Turbine','PV','Grid Purchase','Storage Dis.','Elec. Load','Grid Sale','Storage Ch.','Location','bestoutside');
grid on; xlim([0.5, 24.5]);

% Microgrid 3 Power Balance
subplot(3,2,5);
sources_mg3 = [P_e_GT_mg3', P_e_pv_mg3', P_buy_mg3', P_batd_mg3'];
loads_mg3 = [-L_e_mg3', -P_sell_mg3', -P_batc_mg3'];
bar(time_axis, [sources_mg3, loads_mg3], 'stacked');
title('Microgrid 3 Power Balance');
xlabel('Time (h)');
ylabel('Power (kW)');
legend('Gas Turbine','PV','Grid Purchase','Storage Dis.','Elec. Load','Grid Sale','Storage Ch.','Location','bestoutside');
grid on; xlim([0.5, 24.5]);

fprintf('Visualization finished!\n');
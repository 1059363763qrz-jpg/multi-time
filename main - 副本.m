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
%                  结果可视化代码 (更新版)
%  =====================================================================

final_iter = iter-1; % 保存最终迭代次数
time_axis = 1:T;

% 1. 绘制各主体目标函数值的收敛曲线
figure;
plot(1:final_iter, Obj_DSO(1:final_iter), 'r-o', 'LineWidth', 1.5);
hold on;
plot(1:final_iter, Obj_MG1(1:final_iter), 'g-s', 'LineWidth', 1.5);
plot(1:final_iter, Obj_MG2(1:final_iter), 'b-^', 'LineWidth', 1.5);
plot(1:final_iter, Obj_MG3(1:final_iter), 'm-d', 'LineWidth', 1.5);
plot(1:final_iter, Obj_SESO(1:final_iter), 'k-*', 'LineWidth', 1.5);
xlabel('迭代次数');
ylabel('目标函数值');
title('各主体目标函数值收敛过程');
legend('DSO', 'Microgrid 1', 'Microgrid 2', 'Microgrid 3', 'SESO', 'Location', 'best');
grid on;
box off;

% 2. 绘制残差的收敛曲线
figure;
plot(1:final_iter, Residual(1:final_iter), 'b-x', 'LineWidth', 1.5);
xlabel('迭代次数');
ylabel('残差');
title('ADMM 算法残差收敛过程');
grid on;
box off;

% 3. 绘制DSO与SESO之间的储能交互结果
figure;
bar(time_axis, [P_dso_charge(final_iter+1,:)', -P_dso_discharge(final_iter+1,:)'], 'stacked');
xlabel('时间 (h)');
ylabel('功率 (kW)');
title('DSO 与 SESO 的储能交互功率');
legend('DSO 从 SESO 充电', 'DSO 向 SESO 放电', 'Location', 'northwest');
grid on;
xlim([0.5, 24.5]);

% 4. 绘制各微网与SESO的储能租赁结果
figure;

% 微网1
subplot(3,1,1);
bar(time_axis, [P_mg1_lease_c(final_iter+1,:)', P_mg1_lease_d(final_iter+1,:)']);
ylabel('功率 (kW)');
title('微网1向SESO出租的储能功率');
legend('充电服务', '放电服务');
grid on;
xlim([0.5, 24.5]);

% 微网2
subplot(3,1,2);
bar(time_axis, [P_mg2_lease_c(final_iter+1,:)', P_mg2_lease_d(final_iter+1,:)']);
ylabel('功率 (kW)');
title('微网2向SESO出租的储能功率');
legend('充电服务', '放电服务');
grid on;
xlim([0.5, 24.5]);

% 微网3
subplot(3,1,3);
bar(time_axis, [P_mg3_lease_c(final_iter+1,:)', P_mg3_lease_d(final_iter+1,:)']);
xlabel('时间 (h)');
ylabel('功率 (kW)');
title('微网3向SESO出租的储能功率');
legend('充电服务', '放电服务');
grid on;
xlim([0.5, 24.5]);

% 5. 绘制微网与配电网的净交换功率
figure;
plot(time_axis, P_mg1_net_exchange(final_iter+1,:), 'r-o', 'LineWidth', 1.5);
hold on;
plot(time_axis, P_mg2_net_exchange(final_iter+1,:), 'g-s', 'LineWidth', 1.5);
plot(time_axis, P_mg3_net_exchange(final_iter+1,:), 'b-^', 'LineWidth', 1.5);
yline(0, 'k--', 'LineWidth', 1); % 添加零线
xlabel('时间 (h)');
ylabel('净交换功率 (kW)');
title('各微网与配电网的净交换功率 (正为售电, 负为购电)');
legend('Microgrid 1', 'Microgrid 2', 'Microgrid 3', 'Location', 'best');
grid on;
xlim([1, 24]);
%% =====================================================================
%                  新增：详细优化结果可视化代码 (更新版)
%  =====================================================================

% --- 步骤 1: 使用收敛后的变量，重新运行一次各子问题以获取详细结果 ---
fprintf('正在生成最终可视化结果...\n');

% 获取MG1的详细结果
[~, ~, ~, ~, E_bat_mg1, P_batc_mg1, P_batd_mg1, P_e_GT_mg1, P_e_wd_mg1, P_buy_mg1, P_sell_mg1, L_e_mg1] = Fun_MG1_Lease(...
    P_seso_from_mg1_c(final_iter+1,:), P_seso_from_mg1_d(final_iter+1,:), lambda_mg1_c, lambda_mg1_d, rho);

% 获取MG2的详细结果
[~, ~, ~, ~, E_bat_mg2, P_batc_mg2, P_batd_mg2, P_e_GT_mg2, P_e_pv_mg2, P_buy_mg2, P_sell_mg2, L_e_mg2] = Fun_MG2_Lease(...
    P_seso_from_mg2_c(final_iter+1,:), P_seso_from_mg2_d(final_iter+1,:), lambda_mg2_c, lambda_mg2_d, rho);

% 获取MG3的详细结果
[~, ~, ~, ~, E_bat_mg3, P_batc_mg3, P_batd_mg3, P_e_GT_mg3, P_e_pv_mg3, P_buy_mg3, P_sell_mg3, L_e_mg3] = Fun_MG3_Lease(...
    P_seso_from_mg3_c(final_iter+1,:), P_seso_from_mg3_d(final_iter+1,:), lambda_mg3_c, lambda_mg3_d, rho);

% 获取DSO的详细结果
[~, ~, ~, P_grid_exchange, PLOAD_dso, P_loss_dso, P_pv_dso, P_wt_dso, PLOAD_original, Lshift_effect, SIL_effect] = Fun_DSO(...
    P_seso_to_dso_c(final_iter+1,:), P_seso_to_dso_d(final_iter+1,:), lambda_dso_c, lambda_dso_d, ...
    P_mg1_net_exchange(final_iter+1,:), P_mg2_net_exchange(final_iter+1,:), P_mg3_net_exchange(final_iter+1,:), rho);
P_buy_dso = max(0, P_grid_exchange);
P_sell_dso = -min(0, P_grid_exchange);

% 获取SESO的详细结果
[~, ~, ~, ~, ~, ~, ~, ~, ~, E_seso, P_seso_ch, P_seso_dis] = Fun_SESO(...
    P_dso_charge(final_iter+1,:), P_dso_discharge(final_iter+1,:), lambda_dso_c, lambda_dso_d, ...
    P_mg1_lease_c(final_iter+1,:), P_mg1_lease_d(final_iter+1,:), lambda_mg1_c, lambda_mg1_d, ...
    P_mg2_lease_c(final_iter+1,:), P_mg2_lease_d(final_iter+1,:), lambda_mg2_c, lambda_mg2_d, ...
    P_mg3_lease_c(final_iter+1,:), P_mg3_lease_d(final_iter+1,:), lambda_mg3_c, lambda_mg3_d, rho);

% --- 步骤 2: 绘制储能状态图 ---

% 绘制SESO和各微网的储能状态
figure;
sgtitle('各主体储能系统状态 (收敛结果)');

% SESO 自有储能 SOC
subplot(4,2,1);
plot(time_axis, E_seso(1:T), 'r-s', 'LineWidth', 1.5);
ylabel('电量 (kWh)');
title('SESO 自有储能SOC');
grid on; xlim([1, 24]);

% SESO 自有储能充放电功率
subplot(4,2,2);
bar(time_axis, [P_seso_ch', -P_seso_dis'], 'stacked');
ylabel('功率 (kW)');
title('SESO 自有储能充放电功率');
legend('充电', '放电');
grid on; xlim([0.5, 24.5]);

% 微网1 储能 SOC
subplot(4,2,3);
plot(time_axis, E_bat_mg1, 'g-^', 'LineWidth', 1.5);
ylabel('电量 (kWh)');
title('微网1 储能SOC');
grid on; xlim([1, 24]);

% 微网1 储能充放电功率
subplot(4,2,4);
bar(time_axis, [P_batc_mg1', -P_batd_mg1'], 'stacked');
ylabel('功率 (kW)');
title('微网1 储能充放电功率');
legend('充电', '放电');
grid on; xlim([0.5, 24.5]);

% 微网2 储能 SOC
subplot(4,2,5);
plot(time_axis, E_bat_mg2, 'b-o', 'LineWidth', 1.5);
ylabel('电量 (kWh)');
title('微网2 储能SOC');
grid on; xlim([1, 24]);

% 微网2 储能充放电功率
subplot(4,2,6);
bar(time_axis, [P_batc_mg2', -P_batd_mg2'], 'stacked');
ylabel('功率 (kW)');
title('微网2 储能充放电功率');
legend('充电', '放电');
grid on; xlim([0.5, 24.5]);

% 微网3 储能 SOC
subplot(4,2,7);
plot(time_axis, E_bat_mg3, 'm-d', 'LineWidth', 1.5);
xlabel('时间 (h)');
ylabel('电量 (kWh)');
title('微网3 储能SOC');
grid on; xlim([1, 24]);

% 微网3 储能充放电功率
subplot(4,2,8);
bar(time_axis, [P_batc_mg3', -P_batd_mg3'], 'stacked');
xlabel('时间 (h)');
ylabel('功率 (kW)');
title('微网3 储能充放电功率');
legend('充电', '放电');
grid on; xlim([0.5, 24.5]);

% --- 步骤 3: 绘制功率平衡图 ---
figure;
sgtitle('各主体功率平衡 (源为正，荷为负)');

% DSO 功率平衡
% DSO 功率平衡
subplot(3,2,1);
mg1_source = max(0, P_mg1_net_exchange(final_iter+1,:));
mg1_load = -min(0, P_mg1_net_exchange(final_iter+1,:));
mg2_source = max(0, P_mg2_net_exchange(final_iter+1,:));
mg2_load = -min(0, P_mg2_net_exchange(final_iter+1,:));
mg3_source = max(0, P_mg3_net_exchange(final_iter+1,:));
mg3_load = -min(0, P_mg3_net_exchange(final_iter+1,:));
P_buy_dso = max(0, P_grid_exchange);      % 提取所有正值（购电）
P_sell_dso = -min(0, P_grid_exchange);     % 提取所有负值（售电）并取正
sources_dso = [P_buy_dso', sum(P_pv_dso,1)', sum(P_wt_dso,1)', mg1_source', mg2_source', mg3_source', P_dso_discharge(final_iter+1,:)'];
loads_dso = [-sum(PLOAD_original,1)', sum(Lshift_effect,1)', sum(SIL_effect,1)', ...
             -P_sell_dso', -P_loss_dso', -mg1_load', -mg2_load', -mg3_load', -P_dso_charge(final_iter+1,:)'];

% 绘图
bar(time_axis, [sources_dso, loads_dso], 'stacked');
title('DSO 功率平衡 (计及需求响应分解)');
ylabel('功率 (kW)');

% 更新图例以反映新的负荷分解
legend('主网购电','光伏','风电','从微网购电','SESO放电', ...
       '原始负荷','负荷平移(减少)','负荷削减', ...
       '售电至主网','网损','向微网售电','SESO充电', ...
       'Location','bestoutside');
% ======================== END OF MODIFICATION ========================

grid on; xlim([0.5, 24.5]);


% SESO 功率平衡
subplot(3,2,2);
sources_seso = [P_seso_dis', P_dso_charge(final_iter+1,:)', P_mg1_lease_d(final_iter+1,:)', P_mg2_lease_d(final_iter+1,:)', P_mg3_lease_d(final_iter+1,:)'];
loads_seso = [-P_seso_ch', -P_dso_discharge(final_iter+1,:)', -P_mg1_lease_c(final_iter+1,:)', -P_mg2_lease_c(final_iter+1,:)', -P_mg3_lease_c(final_iter+1,:)'];
bar(time_axis, [sources_seso, loads_seso], 'stacked');
title('SESO 功率平衡');
ylabel('功率 (kW)');
legend('自有放电','DSO充电','MG1租赁放','MG2租赁放','MG3租赁放','自有充电','DSO放电','MG1租赁充','MG2租赁充','MG3租赁充','Location','bestoutside');
grid on; xlim([0.5, 24.5]);

% 微网1 功率平衡
subplot(3,2,3);
sources_mg1 = [P_e_GT_mg1', P_e_wd_mg1', P_buy_mg1', P_batd_mg1'];
loads_mg1 = [-L_e_mg1', -P_sell_mg1', -P_batc_mg1'];
bar(time_axis, [sources_mg1, loads_mg1], 'stacked');
title('微网1 功率平衡');
ylabel('功率 (kW)');
legend('燃气轮机','风电','购电','储能放电','电负荷','售电','储能充电','Location','bestoutside');
grid on; xlim([0.5, 24.5]);

% 微网2 功率平衡
subplot(3,2,4);
sources_mg2 = [P_e_GT_mg2', P_e_pv_mg2', P_buy_mg2', P_batd_mg2'];
loads_mg2 = [-L_e_mg2', -P_sell_mg2', -P_batc_mg2'];
bar(time_axis, [sources_mg2, loads_mg2], 'stacked');
title('微网2 功率平衡');
ylabel('功率 (kW)');
legend('燃气轮机','光伏','购电','储能放电','电负荷','售电','储能充电','Location','bestoutside');
grid on; xlim([0.5, 24.5]);

% 微网3 功率平衡
subplot(3,2,5);
sources_mg3 = [P_e_GT_mg3', P_e_pv_mg3', P_buy_mg3', P_batd_mg3'];
loads_mg3 = [-L_e_mg3', -P_sell_mg3', -P_batc_mg3'];
bar(time_axis, [sources_mg3, loads_mg3], 'stacked');
title('微网3 功率平衡');
xlabel('时间 (h)');
ylabel('功率 (kW)');
legend('燃气轮机','光伏','购电','储能放电','电负荷','售电','储能充电','Location','bestoutside');
grid on; xlim([0.5, 24.5]);

fprintf('可视化结果生成完毕！\n');
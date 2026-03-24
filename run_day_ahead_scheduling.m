function [plan_DA, admm_results_DA] = run_day_ahead_scheduling(dso_data_da, mg_data_da)
    % 执行日前调度
    % V2.0: 增加返回 MG GT 的日前计划

    %% ADMM 迭代参数
    maxIter = 200; rho = 1; tolerant = 2; iter = 1;
    T = 24;
    admm_results_DA = struct();

    % === 初始化 Primal 变量 ===
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

    % === 初始化 Dual 变量 ===
    lambda_dso_c = zeros(1, T); lambda_dso_d = zeros(1, T);
    lambda_mg1_c = zeros(1, T); lambda_mg1_d = zeros(1, T);
    lambda_mg2_c = zeros(1, T); lambda_mg2_d = zeros(1, T);
    lambda_mg3_c = zeros(1, T); lambda_mg3_d = zeros(1, T);

    % === 历史数据记录 ===
    Obj_DSO = []; Obj_MG1 = []; Obj_MG2 = []; Obj_MG3 = []; Obj_SESO = [];
    Residual = [];

    %% ADMM 迭代主循环
    while iter <= maxIter
        fprintf('正在进行第 %d 次迭代...  %f\n', iter);
        % --- 步骤 1: 求解微网子问题 ---
        [P_mg1_net_exchange(iter+1,:), P_mg1_lease_c(iter+1,:), P_mg1_lease_d(iter+1,:), obj_mg1] = Fun_MG1_Lease(...
            P_seso_from_mg1_c(iter,:), P_seso_from_mg1_d(iter,:), lambda_mg1_c, lambda_mg1_d, rho);
        Obj_MG1(iter) = obj_mg1;

        [P_mg2_net_exchange(iter+1,:), P_mg2_lease_c(iter+1,:), P_mg2_lease_d(iter+1,:), obj_mg2] = Fun_MG2_Lease(...
            P_seso_from_mg2_c(iter,:), P_seso_from_mg2_d(iter,:), lambda_mg2_c, lambda_mg2_d, rho);
        Obj_MG2(iter) = obj_mg2;

        [P_mg3_net_exchange(iter+1,:), P_mg3_lease_c(iter+1,:), P_mg3_lease_d(iter+1,:), obj_mg3] = Fun_MG3_Lease(...
            P_seso_from_mg3_c(iter,:), P_seso_from_mg3_d(iter,:), lambda_mg3_c, lambda_mg3_d, rho);
        Obj_MG3(iter) = obj_mg3;

        % --- 步骤 2: 求解DSO子问题 ---
        [P_dso_charge(iter+1,:), P_dso_discharge(iter+1,:), obj_dso] = Fun_DSO(...
            P_seso_to_dso_c(iter,:), P_seso_to_dso_d(iter,:), lambda_dso_c, lambda_dso_d, ...
            P_mg1_net_exchange(iter+1,:), P_mg2_net_exchange(iter+1,:), P_mg3_net_exchange(iter+1,:), rho);
        Obj_DSO(iter) = obj_dso;

        % --- 步骤 3: 求解SESO子问题 ---
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

        % --- 步骤 4: 更新拉格朗日乘子 ---
        lambda_dso_c = lambda_dso_c + rho * (P_dso_charge(iter+1,:) - P_seso_to_dso_c(iter+1,:));
        lambda_dso_d = lambda_dso_d + rho * (P_dso_discharge(iter+1,:) - P_seso_to_dso_d(iter+1,:));
        lambda_mg1_c = lambda_mg1_c + rho * (P_mg1_lease_c(iter+1,:) - P_seso_from_mg1_c(iter+1,:));
        lambda_mg1_d = lambda_mg1_d + rho * (P_mg1_lease_d(iter+1,:) - P_seso_from_mg1_d(iter+1,:));
        lambda_mg2_c = lambda_mg2_c + rho * (P_mg2_lease_c(iter+1,:) - P_seso_from_mg2_c(iter+1,:));
        lambda_mg2_d = lambda_mg2_d + rho * (P_mg2_lease_d(iter+1,:) - P_seso_from_mg2_d(iter+1,:));
        lambda_mg3_c = lambda_mg3_c + rho * (P_mg3_lease_c(iter+1,:) - P_seso_from_mg3_c(iter+1,:));
        lambda_mg3_d = lambda_mg3_d + rho * (P_mg3_lease_d(iter+1,:) - P_seso_from_mg3_d(iter+1,:));

        % --- 步骤 5: 计算并检查收敛性 ---
        res_dso_seso = norm(P_dso_charge(iter+1,:) - P_seso_to_dso_c(iter+1,:)) + norm(P_dso_discharge(iter+1,:) - P_seso_to_dso_d(iter+1,:));
        res_mg1_seso = norm(P_mg1_lease_c(iter+1,:) - P_seso_from_mg1_c(iter+1,:)) + norm(P_mg1_lease_d(iter+1,:) - P_seso_from_mg1_d(iter+1,:));
        res_mg2_seso = norm(P_mg2_lease_c(iter+1,:) - P_seso_from_mg2_c(iter+1,:)) + norm(P_mg2_lease_d(iter+1,:) - P_seso_from_mg2_d(iter+1,:));
        res_mg3_seso = norm(P_mg3_lease_c(iter+1,:) - P_seso_from_mg3_c(iter+1,:)) + norm(P_mg3_lease_d(iter+1,:) - P_seso_from_mg3_d(iter+1,:));

        current_residual = res_dso_seso + res_mg1_seso + res_mg2_seso + res_mg3_seso;
        Residual(iter) = current_residual;
        fprintf('当前残差: %f\n', current_residual);
        if iter == maxIter
            disp('已达到最大迭代次数，算法未收敛。');
        end

        iter = iter + 1;

        if current_residual < tolerant
            fprintf('求解成功，正在整理最终结果...\n');
            
            % --- 获取详细的最终结果 ---
            % ======================== MODIFICATION START ========================
            % 获取MG1的详细结果 (增加捕获 P_e_GT_val)
            [~, ~, ~, ~, E_bat_val1, P_batc_val1, P_batd_val1, P_e_GT_val1, ~, ~, ~, ~, P_e_cut_val1, P_e_tran_val1, ~] = Fun_MG1_Lease(...
                P_seso_from_mg1_c(iter,:), P_seso_from_mg1_d(iter,:), lambda_mg1_c, lambda_mg1_d, rho);

            % 获取MG2的详细结果 (增加捕获 P_e_GT_val)
            [~, ~, ~, ~, E_bat_val2, P_batc_val2, P_batd_val2, P_e_GT_val2, ~, ~, ~, ~, P_e_cut_val2, P_e_tran_val2, ~] = Fun_MG2_Lease(...
                P_seso_from_mg2_c(iter,:), P_seso_from_mg2_d(iter,:), lambda_mg2_c, lambda_mg2_d, rho);

            % 获取MG3的详细结果 (增加捕获 P_e_GT_val)
            [~, ~, ~, ~, E_bat_val3, P_batc_val3, P_batd_val3, P_e_GT_val3, ~, ~, ~, ~, P_e_cut_val3, P_e_tran_val3, ~] = Fun_MG3_Lease(...
                P_seso_from_mg3_c(iter,:), P_seso_from_mg3_d(iter,:), lambda_mg3_c, lambda_mg3_d, rho);
            % ======================== MODIFICATION END ========================

            % 获取DSO的详细结果
            [~, ~, ~, P_grid_exchange, ~, P_loss_dso, ~, ~, ~, Lshift_effect, SIL_effect] = Fun_DSO(...
                P_seso_to_dso_c(iter,:), P_seso_to_dso_d(iter,:), lambda_dso_c, lambda_dso_d, ...
                P_mg1_net_exchange(iter,:), P_mg2_net_exchange(iter,:), P_mg3_net_exchange(iter,:), rho);
            
            % 获取SESO的详细结果
            [~, ~, ~, ~, ~, ~, ~, ~, ~, E_seso, P_seso_ch, P_seso_dis] = Fun_SESO(...
                P_dso_charge(iter,:), P_dso_discharge(iter,:), lambda_dso_c, lambda_dso_d, ...
                P_mg1_lease_c(iter,:), P_mg1_lease_d(iter,:), lambda_mg1_c, lambda_mg1_d, ...
                P_mg2_lease_c(iter,:), P_mg2_lease_d(iter,:), lambda_mg2_c, lambda_mg2_d, ...
                P_mg3_lease_c(iter,:), P_mg3_lease_d(iter,:), lambda_mg3_c, lambda_mg3_d, rho);

            % --- 赋值给 plan_DA ---
            plan_DA.P_dso_grid = P_grid_exchange;
            plan_DA.P_dso_charge = P_dso_charge(iter-1,:);
            plan_DA.P_dso_discharge = P_dso_discharge(iter-1,:);
            plan_DA.P_loss = P_loss_dso;
            plan_DA.Lshift_effect_dso = Lshift_effect;
            plan_DA.SIL_effect_dso = SIL_effect;
            
            plan_DA.P_mg1_net_exchange = P_mg1_net_exchange(iter-1,:);
            plan_DA.P_mg1_lease_c = P_mg1_lease_c(iter-1,:);
            plan_DA.P_mg1_lease_d = P_mg1_lease_d(iter-1,:);
            plan_DA.P_e_cut_mg1 = P_e_cut_val1;
            plan_DA.P_e_tran_mg1 = P_e_tran_val1;
            plan_DA.E_mg1 = E_bat_val1;
            plan_DA.P_mg1_batc = P_batc_val1;
            plan_DA.P_mg1_batd = P_batd_val1;
            % ======================== MODIFICATION START ========================
            plan_DA.P_mg1_gt = P_e_GT_val1; % 新增
            % ======================== MODIFICATION END ========================
            
            plan_DA.P_mg2_net_exchange = P_mg2_net_exchange(iter-1,:);
            plan_DA.P_mg2_lease_c = P_mg2_lease_c(iter-1,:);
            plan_DA.P_mg2_lease_d = P_mg2_lease_d(iter-1,:);
            plan_DA.P_e_cut_mg2 = P_e_cut_val2;
            plan_DA.P_e_tran_mg2 = P_e_tran_val2;
            plan_DA.E_mg2 = E_bat_val2;
            plan_DA.P_mg2_batc = P_batc_val2;
            plan_DA.P_mg2_batd = P_batd_val2;
            % ======================== MODIFICATION START ========================
            plan_DA.P_mg2_gt = P_e_GT_val2; % 新增
            % ======================== MODIFICATION END ========================

            plan_DA.P_mg3_net_exchange = P_mg3_net_exchange(iter-1,:);
            plan_DA.P_mg3_lease_c = P_mg3_lease_c(iter-1,:);
            plan_DA.P_mg3_lease_d = P_mg3_lease_d(iter-1,:);
            plan_DA.P_e_cut_mg3 = P_e_cut_val3;
            plan_DA.P_e_tran_mg3 = P_e_tran_val3;
            plan_DA.E_mg3 = E_bat_val3;
            plan_DA.P_mg3_batc = P_batc_val3;
            plan_DA.P_mg3_batd = P_batd_val3;
            % ======================== MODIFICATION START ========================
            plan_DA.P_mg3_gt = P_e_GT_val3; % 新增
            % ======================== MODIFICATION END ========================

            plan_DA.E_seso = E_seso;
            plan_DA.P_seso_ch = P_seso_ch;
            plan_DA.P_seso_dis = P_seso_dis;
            
            % --- 存储ADMM结果 ---
            admm_results_DA.Obj_DSO = Obj_DSO(1:iter-1);
            admm_results_DA.Obj_MG1 = Obj_MG1(1:iter-1);
            admm_results_DA.Obj_MG2 = Obj_MG2(1:iter-1);
            admm_results_DA.Obj_MG3 = Obj_MG3(1:iter-1);
            admm_results_DA.Obj_SESO = Obj_SESO(1:iter-1);
            admm_results_DA.Residual = Residual(1:iter-1);
            break; % 收敛，跳出循环
        else
            % 如果未达到收敛容差且达到最大迭代次数
            if iter > maxIter
                 plan_DA = []; % 返回空结果表示失败
                 admm_results_DA = [];
            end
        end % end if current_residual < tolerant
    end % end while
end
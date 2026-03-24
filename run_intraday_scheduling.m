function [deltas, success] = run_intraday_scheduling(plan_DA, forecasts_RT, k, H)
    % 执行日内滚动调度 ADMM 算法 (基于调整量)
    % 版本更新: 修正了最终结果的索引错误

    %% ADMM 参数
    maxIter = 100; 
    rho = 0.1;       % 罚因子 (日内可以适当增大以加速收敛, 原10可能过大，调整为2)
    tolerant = 5;   % 收敛精度
    iter = 1;

    % =================================================================
    %      初始化 Primal 和 Dual 变量 (全部为调整量 Delta)
    % =================================================================
    
    % --- Primal 变量 (迭代中的调整量) ---
    delta_dso_seso_c = zeros(maxIter+1, H); delta_dso_seso_d = zeros(maxIter+1, H);
    delta_seso_dso_c = zeros(maxIter+1, H); delta_seso_dso_d = zeros(maxIter+1, H);
    delta_mg1_seso_c = zeros(maxIter+1, H); delta_mg1_seso_d = zeros(maxIter+1, H);
    delta_seso_mg1_c = zeros(maxIter+1, H); delta_seso_mg1_d = zeros(maxIter+1, H);
    delta_mg2_seso_c = zeros(maxIter+1, H); delta_mg2_seso_d = zeros(maxIter+1, H);
    delta_seso_mg2_c = zeros(maxIter+1, H); delta_seso_mg2_d = zeros(maxIter+1, H);
    delta_mg3_seso_c = zeros(maxIter+1, H); delta_mg3_seso_d = zeros(maxIter+1, H);
    delta_seso_mg3_c = zeros(maxIter+1, H); delta_seso_mg3_d = zeros(maxIter+1, H);
    delta_mg1_net = zeros(maxIter+1, H);
    delta_mg2_net = zeros(maxIter+1, H);
    delta_mg3_net = zeros(maxIter+1, H);
    delta_dso_grid = zeros(maxIter+1, H); % <<<<<<< 新增：需要初始化dso_grid的存储
    
    % --- Dual 变量 (拉格朗日乘子) ---
    lambda_dso_c = zeros(1, H); lambda_dso_d = zeros(1, H);
    lambda_mg1_c = zeros(1, H); lambda_mg1_d = zeros(1, H);
    lambda_mg2_c = zeros(1, H); lambda_mg2_d = zeros(1, H);
    lambda_mg3_c = zeros(1, H); lambda_mg3_d = zeros(1, H);
    
    %% 日内 ADMM 循环
    while iter <= maxIter
        fprintf('--- 日内 ADMM 迭代: %d ---\n', iter);
        
        % 步骤 1: 求解微网子问题
        [delta_mg1_net(iter+1,:), delta_mg1_seso_c(iter+1,:), delta_mg1_seso_d(iter+1,:)] = Fun_MG1_intraday(plan_DA, forecasts_RT.mg, k, H, ...
            delta_seso_mg1_c(iter,:), delta_seso_mg1_d(iter,:), lambda_mg1_c, lambda_mg1_d, rho);
        
        [delta_mg2_net(iter+1,:), delta_mg2_seso_c(iter+1,:), delta_mg2_seso_d(iter+1,:)] = Fun_MG2_intraday(plan_DA, forecasts_RT.mg, k, H, ...
            delta_seso_mg2_c(iter,:), delta_seso_mg2_d(iter,:), lambda_mg2_c, lambda_mg2_d, rho);
        
        [delta_mg3_net(iter+1,:), delta_mg3_seso_c(iter+1,:), delta_mg3_seso_d(iter+1,:)] = Fun_MG3_intraday(plan_DA, forecasts_RT.mg, k, H, ...
            delta_seso_mg3_c(iter,:), delta_seso_mg3_d(iter,:), lambda_mg3_c, lambda_mg3_d, rho);

        % 步骤 2: 求解DSO子问题
        [delta_dso_seso_c(iter+1,:), delta_dso_seso_d(iter+1,:), delta_dso_grid(iter+1,:)] = Fun_DSO_intraday(plan_DA, forecasts_RT.dso, k, H, ...
            delta_seso_dso_c(iter,:), delta_seso_dso_d(iter,:), lambda_dso_c, lambda_dso_d, ...
            delta_mg1_net(iter+1,:), delta_mg2_net(iter+1,:), delta_mg3_net(iter+1,:), rho);

        % 步骤 3: 求解SESO子问题
        [delta_seso_dso_c(iter+1,:), delta_seso_dso_d(iter+1,:), ...
         delta_seso_mg1_c(iter+1,:), delta_seso_mg1_d(iter+1,:), ...
         delta_seso_mg2_c(iter+1,:), delta_seso_mg2_d(iter+1,:), ...
         delta_seso_mg3_c(iter+1,:), delta_seso_mg3_d(iter+1,:)] = Fun_SESO_intraday(plan_DA, k, H, ...
            delta_dso_seso_c(iter+1,:), delta_dso_seso_d(iter+1,:), lambda_dso_c, lambda_dso_d, ...
            delta_mg1_seso_c(iter+1,:), delta_mg1_seso_d(iter+1,:), lambda_mg1_c, lambda_mg1_d, ...
            delta_mg2_seso_c(iter+1,:), delta_mg2_seso_d(iter+1,:), lambda_mg2_c, lambda_mg2_d, ...
            delta_mg3_seso_c(iter+1,:), delta_mg3_seso_d(iter+1,:), lambda_mg3_c, lambda_mg3_d, rho);

        % 步骤 4: 更新拉格朗日乘子
        lambda_dso_c = lambda_dso_c + rho * (delta_dso_seso_c(iter+1,:) - delta_seso_dso_c(iter+1,:));
        lambda_dso_d = lambda_dso_d + rho * (delta_dso_seso_d(iter+1,:) - delta_seso_dso_d(iter+1,:));
        lambda_mg1_c = lambda_mg1_c + rho * (delta_mg1_seso_c(iter+1,:) - delta_seso_mg1_c(iter+1,:));
        lambda_mg1_d = lambda_mg1_d + rho * (delta_mg1_seso_d(iter+1,:) - delta_seso_mg1_d(iter+1,:));
        lambda_mg2_c = lambda_mg2_c + rho * (delta_mg2_seso_c(iter+1,:) - delta_seso_mg2_c(iter+1,:));
        lambda_mg2_d = lambda_mg2_d + rho * (delta_mg2_seso_d(iter+1,:) - delta_seso_mg2_d(iter+1,:));
        lambda_mg3_c = lambda_mg3_c + rho * (delta_mg3_seso_c(iter+1,:) - delta_seso_mg3_c(iter+1,:));
        lambda_mg3_d = lambda_mg3_d + rho * (delta_mg3_seso_d(iter+1,:) - delta_seso_mg3_d(iter+1,:));

        % 步骤 5: 检查收敛性
        res = norm(delta_dso_seso_c(iter+1,:) - delta_seso_dso_c(iter+1,:)) + ...
              norm(delta_dso_seso_d(iter+1,:) - delta_seso_dso_d(iter+1,:)) + ...
              norm(delta_mg1_seso_c(iter+1,:) - delta_seso_mg1_c(iter+1,:)) + ...
              norm(delta_mg1_seso_d(iter+1,:) - delta_seso_mg1_d(iter+1,:)) + ...
              norm(delta_mg2_seso_c(iter+1,:) - delta_seso_mg2_c(iter+1,:)) + ...
              norm(delta_mg2_seso_d(iter+1,:) - delta_seso_mg2_d(iter+1,:)) + ...
              norm(delta_mg3_seso_c(iter+1,:) - delta_seso_mg3_c(iter+1,:)) + ...
              norm(delta_mg3_seso_d(iter+1,:) - delta_seso_mg3_d(iter+1,:));

        fprintf('    当前残差: %f\n', res);

        if res < tolerant
            fprintf('日内调度在第 %d 次迭代收敛！\n', iter);
            break; 
        end
        if iter == maxIter, disp('日内调度已达到最大迭代次数，未收敛。'); end
        iter = iter + 1;
    end
    
    % --- 输出最终的调整量 ---
    if res < tolerant
        success = true;
        final_iter = iter; % 使用当前iter作为最终迭代次数
        % =================== 关键修正：所有索引统一为 final_iter+1 ===================
        deltas.dso_seso_c = delta_dso_seso_c(final_iter+1,:);
        deltas.dso_seso_d = delta_dso_seso_d(final_iter+1,:);
        deltas.dso_grid   = delta_dso_grid(final_iter+1,:);
        deltas.mg1_net    = delta_mg1_net(final_iter+1,:);
        deltas.mg2_net    = delta_mg2_net(final_iter+1,:);
        deltas.mg3_net    = delta_mg3_net(final_iter+1,:);
        deltas.mg1_seso_c = delta_mg1_seso_c(final_iter+1,:);
        deltas.mg1_seso_d = delta_mg1_seso_d(final_iter+1,:);
        deltas.mg2_seso_c = delta_mg2_seso_c(final_iter+1,:);
        deltas.mg2_seso_d = delta_mg2_seso_d(final_iter+1,:);
        deltas.mg3_seso_c = delta_mg3_seso_c(final_iter+1,:);
        deltas.mg3_seso_d = delta_mg3_seso_d(final_iter+1,:);
        % =======================================================================
    else
        success = false;
        deltas = [];
    end
end
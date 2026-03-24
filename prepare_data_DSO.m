function dso_data = prepare_data_DSO(T_DA, T_RT)
    % 为DSO准备日前预测和日内实际数据
    
   mpc = IEEE33BW;
    PLOAD_DA = mpc.Pload * 100; % 节点有功负荷 (kW)
    T = 24; nwt = 2; npv = 1;

    % 生成可再生能源预测
    v=xlsread('风速数据.xlsx','B2:B25');
    v_in=3; v_r=11.3; v_out=25; S_wt=[1; 1.5];
    P_wt_forecast_DA = zeros(nwt, T);
    for i=1:1:nwt
        for t=1:1:T
          if (v(t) >= v_in) && (v(t) < v_r)
              P_wt_forecast_DA(i,t)=[(v(t)-v_in)/(v_r-v_in)]*S_wt(i);
          elseif (v(t) >= v_r) && (v(t) < v_out)
              P_wt_forecast_DA(i,t)=S_wt(i);
          end
        end
    end

    RR=xlsread('光照强度数据.xlsx','B2:B25');
    R_STD=1000; R_C=150; S_pv=1.5;
    P_pv_forecast_DA = zeros(npv, T);
    for i=1:1:npv
        for t=1:1:T
            if(RR(t) > 0) && (RR(t) <= R_C)
                P_pv_forecast_DA(i,t)=S_pv*((RR(t))^2)/(R_STD*R_C);
            elseif (RR(t) > R_C) && (RR(t) <= R_STD)
                P_pv_forecast_DA(i,t)=S_pv*RR(t)/R_STD;
            elseif (RR(t) > R_STD)
                P_pv_forecast_DA(i,t)=S_pv;
            end
        end
    end

    % 2. 生成日内“真实”数据 (96点, 15分钟/点)
    

    % 对矩阵（负荷、风电）进行插值，需要先转置，插值后再转置回来
    PLOAD_RT = interp1(T_DA, PLOAD_DA', T_RT, 'spline')';
    P_wt_RT   = interp1(T_DA, P_wt_forecast_DA', T_RT, 'spline')';
    
    % --- Bug修复 ---
    % 对于向量（如光伏数据），直接进行插值，无需转置
    P_pv_RT   = interp1(T_DA, P_pv_forecast_DA, T_RT, 'spline');
    % --- 结束修复 ---
    
    % 增加随机噪声
    PLOAD_RT = PLOAD_RT + 0.005 * mean(PLOAD_DA, 'all') * randn(size(PLOAD_RT));
    P_wt_RT = P_wt_RT + 0.0001 * mean(P_wt_forecast_DA, 'all') * randn(size(P_wt_RT));
    P_pv_RT = P_pv_RT + 0.0001 * mean(P_pv_forecast_DA, 'all') * randn(size(P_pv_RT));

    % 确保非负
    PLOAD_RT(PLOAD_RT < 0) = 0;
    P_wt_RT(P_wt_RT < 0) = 0;
    P_pv_RT(P_pv_RT < 0) = 0;

    % 3. 将数据打包成结构体输出
    forecasts_DA.PLOAD = PLOAD_DA;
    forecasts_DA.P_wt  = P_wt_forecast_DA;
    forecasts_DA.P_pv  = P_pv_forecast_DA;

    forecasts_RT.PLOAD = PLOAD_RT;
    forecasts_RT.P_wt  = P_wt_RT;
    forecasts_RT.P_pv  = P_pv_RT;

    dso_data.forecasts_DA = forecasts_DA;
    dso_data.forecasts_RT = forecasts_RT;
end
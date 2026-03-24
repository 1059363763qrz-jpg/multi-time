function mg_data = prepare_data_MGs(T_DA, T_RT)

    % 为微网准备日前预测和日内实际数据
   
    % --- 日前数据 (小时级) ---
    % MG1
    mg_data.forecasts_DA.L_e0_mg1 = [3482,3513,3405,3529,3591,3482,3622,4504,5294,7074,7739,8328,8900,8993,9040,8467,7801,7120,6176,5325,4334,3452,3436,3405]/10;
    mg_data.forecasts_DA.L_h0_mg1 = [1666,1643,1697,1890,2060,2515,4290,4652,4683,4891,4722,4691,4583,4699,4583,4621,4220,4328,4429,4212,3711,2237,1851,1697]/10;
    mg_data.forecasts_DA.Predict_wt_mg1 = [4610,4423,3987,3901,3816,3481,3325,3146,3255,3325,3380,3271,3317,3200,3208,3442,3263,3107,3208,3496,3629,3652,4376,4595]/10;
    % MG2 (示例数据)
    mg_data.forecasts_DA.L_e0_mg2 = [1774,1450,1296,1219,1095,1265,1481,1944,2484,2083,1651,1188,1080,1126,1033,1033,941,1450,2283,3148,3904,3719,2746,2453]/10;
    mg_data.forecasts_DA.L_h0_mg2 = [1666,1643,1697,1890,2060,2515,4290,4652,4683,4891,4722,4691,4583,4699,4583,4621,4220,4328,4429,4212,3711,2237,1851,1697]/10;
    mg_data.forecasts_DA.Predict_pv_mg2 = [0,0,0,0,0,0,967,1287,1583,1833,1918,1942,2004,1957,1669,1076,655,0,0,0,0,0,0,0]/10;
    % MG3 (示例数据)
    mg_data.forecasts_DA.L_e0_mg3 = [1180,1073,1196,1165,1165,1165,1196,1503,1733,2147,2193,2407,2469,2653,2852,2208,2285,2300,3220,3220,2423,1993,1779,1518]/10;
    mg_data.forecasts_DA.L_h0_mg3 = [1494,1448,1325,1302,1317,1494,1594,1833,2080,2211,2311,2388,2303,2526,2434,2326,2164,2126,2118,2519,1841,1625,1494,1394]/10;
    mg_data.forecasts_DA.Predict_pv_mg3 = [0,0,0,0,0,0,663,1084,1903,2277,2386,2480,2402,2168,2012,1474,998,0,0,0,0,0,0,0]/10;

    % --- 日内数据 (15分钟级, 模拟生成) ---
    P_e_load_mg1_base_RT = interp1(T_DA, mg_data.forecasts_DA.L_e0_mg1, T_RT, 'spline');
    P_h_load_mg1_base_RT = interp1(T_DA, mg_data.forecasts_DA.L_h0_mg1, T_RT, 'spline');
    P_wt_mg1_base_RT   = interp1(T_DA, mg_data.forecasts_DA.Predict_wt_mg1, T_RT, 'spline');
    P_wt_mg1_base_RT(P_wt_mg1_base_RT < 0) = 0;
    mg1_e_load_noise = 0.005 * mean(mg_data.forecasts_DA.L_e0_mg1) * randn(1, 96);
    mg1_h_load_noise = 0.005 * mean(mg_data.forecasts_DA.L_h0_mg1) * randn(1, 96);
    mg1_wt_noise   = 0.005 * mean(mg_data.forecasts_DA.Predict_wt_mg1) * randn(1, 96);
    P_e_load_mg1_actual_RT = P_e_load_mg1_base_RT + mg1_e_load_noise;
    P_h_load_mg1_actual_RT = P_h_load_mg1_base_RT + mg1_h_load_noise;
    P_wt_mg1_actual_RT   = P_wt_mg1_base_RT + mg1_wt_noise;
    P_wt_mg1_actual_RT(P_wt_mg1_actual_RT < 0) = 0;
    mg_data.forecasts_RT.L_e0_mg1 = P_e_load_mg1_actual_RT;
    mg_data.forecasts_RT.L_h0_mg1 = P_h_load_mg1_actual_RT;
    mg_data.forecasts_RT.Predict_wt_mg1 = P_wt_mg1_actual_RT;
    
    P_e_load_mg2_base_RT = interp1(T_DA, mg_data.forecasts_DA.L_e0_mg2, T_RT, 'spline');
    P_h_load_mg2_base_RT = interp1(T_DA, mg_data.forecasts_DA.L_h0_mg2, T_RT, 'spline');
    P_pv_mg2_base_RT   = interp1(T_DA, mg_data.forecasts_DA.Predict_pv_mg2, T_RT, 'spline');
    P_pv_mg2_base_RT(P_pv_mg2_base_RT < 0) = 0;
    mg2_e_load_noise = 0.005 * mean(mg_data.forecasts_DA.L_e0_mg2) * randn(1, 96);
    mg2_h_load_noise = 0.005 * mean(mg_data.forecasts_DA.L_h0_mg2) * randn(1, 96);
    mg2_pv_noise   = 0.005 * mean(mg_data.forecasts_DA.Predict_pv_mg2) * randn(1, 96);
    P_e_load_mg2_actual_RT = P_e_load_mg2_base_RT + mg2_e_load_noise;
    P_h_load_mg2_actual_RT = P_h_load_mg2_base_RT + mg2_h_load_noise;
    P_pv_mg2_actual_RT   = P_pv_mg2_base_RT + mg2_pv_noise;
    P_pv_mg2_actual_RT(P_pv_mg2_actual_RT < 0) = 0;
    mg_data.forecasts_RT.L_e0_mg2 = P_e_load_mg2_actual_RT;
    mg_data.forecasts_RT.L_h0_mg2 = P_h_load_mg2_actual_RT;
    mg_data.forecasts_RT.Predict_pv_mg2 = P_pv_mg2_actual_RT;

    P_e_load_mg3_base_RT = interp1(T_DA, mg_data.forecasts_DA.L_e0_mg3, T_RT, 'spline');
    P_h_load_mg3_base_RT = interp1(T_DA, mg_data.forecasts_DA.L_h0_mg3, T_RT, 'spline');
    P_pv_mg3_base_RT   = interp1(T_DA, mg_data.forecasts_DA.Predict_pv_mg3, T_RT, 'spline');
    P_pv_mg3_base_RT(P_pv_mg3_base_RT < 0) = 0;
    mg3_e_load_noise = 0.005 * mean(mg_data.forecasts_DA.L_e0_mg3) * randn(1, 96);
    mg3_h_load_noise = 0.005 * mean(mg_data.forecasts_DA.L_h0_mg3) * randn(1, 96);
    mg3_pv_noise   = 0.005 * mean(mg_data.forecasts_DA.Predict_pv_mg3) * randn(1, 96);
    P_e_load_mg3_actual_RT = P_e_load_mg3_base_RT + mg3_e_load_noise;
    P_h_load_mg3_actual_RT = P_h_load_mg3_base_RT + mg3_h_load_noise;
    P_pv_mg3_actual_RT   = P_pv_mg3_base_RT + mg3_pv_noise;
    P_pv_mg3_actual_RT(P_pv_mg3_actual_RT < 0) = 0;
    mg_data.forecasts_RT.L_e0_mg3 = P_e_load_mg3_actual_RT;
    mg_data.forecasts_RT.L_h0_mg3 = P_h_load_mg3_actual_RT;
    mg_data.forecasts_RT.Predict_pv_mg3 = P_pv_mg3_actual_RT;
    
end
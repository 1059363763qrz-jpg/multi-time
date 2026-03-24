function zero_deltas = get_zero_deltas(window_len)
% get_zero_deltas: 生成所有调整量为0的结构体。

    zero_deltas.dso_grid = zeros(1, window_len);
    zero_deltas.dso_seso_c = zeros(1, window_len);
    zero_deltas.dso_seso_d = zeros(1, window_len);
    zero_deltas.mg1_net = zeros(1, window_len);
    zero_deltas.mg2_net = zeros(1, window_len);
    zero_deltas.mg3_net = zeros(1, window_len);
    zero_deltas.mg1_seso_c = zeros(1, window_len);
    zero_deltas.mg1_seso_d = zeros(1, window_len);
    zero_deltas.mg2_seso_c = zeros(1, window_len);
    zero_deltas.mg2_seso_d = zeros(1, window_len);
    zero_deltas.mg3_seso_c = zeros(1, window_len);
    zero_deltas.mg3_seso_d = zeros(1, window_len);

    % 优化失败时，储能状态由上层逻辑处理
    zero_deltas.E_seso_actual_t1 = NaN;
    zero_deltas.P_seso_ch_actual_t1 = 0;
    zero_deltas.P_seso_dis_actual_t1 = 0;
    zero_deltas.E_mg1_actual_t1 = NaN;
    zero_deltas.P_mg1_batc_actual_t1 = 0;
    zero_deltas.P_mg1_batd_actual_t1 = 0;
    zero_deltas.E_mg2_actual_t1 = NaN;
    zero_deltas.P_mg2_batc_actual_t1 = 0;
    zero_deltas.P_mg2_batd_actual_t1 = 0;
    zero_deltas.E_mg3_actual_t1 = NaN;
    zero_deltas.P_mg3_batc_actual_t1 = 0;
    zero_deltas.P_mg3_batd_actual_t1 = 0;
end

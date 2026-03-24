function scenario_file = generate_fixed_rt_scenario(varargin)
% generate_fixed_rt_scenario: 生成固定RT场景并保存。
% 用法:
%   generate_fixed_rt_scenario();
%   generate_fixed_rt_scenario('rt_scenario_fixed.mat', 2026, T_DA, T_RT);

    scenario_file = 'rt_scenario_fixed.mat';
    rng_seed = 2026;
    T_DA = 1:24;
    T_RT = (15/60):(15/60):24;

    if nargin >= 1 && ~isempty(varargin{1}), scenario_file = varargin{1}; end
    if nargin >= 2 && ~isempty(varargin{2}), rng_seed = varargin{2}; end
    if nargin >= 3 && ~isempty(varargin{3}), T_DA = varargin{3}; end
    if nargin >= 4 && ~isempty(varargin{4}), T_RT = varargin{4}; end

    rng(rng_seed);
    dso_data = prepare_data_DSO(T_DA, T_RT);
    mg_data = prepare_data_MGs(T_DA, T_RT);

    save(scenario_file, 'dso_data', 'mg_data', 'T_DA', 'T_RT', 'rng_seed');
    fprintf('固定日内场景已保存到 %s (rng=%d)。\n', scenario_file, rng_seed);
end

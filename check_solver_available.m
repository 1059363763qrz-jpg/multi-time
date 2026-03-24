function check_solver_available()
% check_solver_available: 前置检查 Gurobi/YALMIP 可用性，失败即停。

    if exist('optimize', 'file') ~= 2 && exist('solvesdp', 'file') ~= 2
        error('未检测到 YALMIP (optimize/solvesdp)。请先将 YALMIP 加入 MATLAB path。');
    end

    has_gurobi = (exist('gurobi', 'file') == 2) || ~isempty(which('gurobi'));
    if ~has_gurobi
        error('Solver not found: gurobi。请先配置 Gurobi 并加入 MATLAB path。');
    end
end

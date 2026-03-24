
function [P_mg1_net_exchange, P_mg1_lease_charge, P_mg1_lease_discharge, obj_mg1, E_bat_val, P_batc_val, P_batd_val, P_e_GT_val, P_e_wd_val, P_buy_val, P_sell_val, L_e_val,P_e_cut_val,P_e_tran_val,P_h_DR_val]  = Fun_MG1_Lease( P_seso_from_mg1_c, P_seso_from_mg1_d, lambda_mg1_c, lambda_mg1_d, rho)
% % Fun_MG_Lease: 考虑储能租赁的微网模型
%微网参数设置,
 % clc; clear; close all;
L_e0=[3482,3513,3405,3529,3591,3482,3622,4504,5294,7074,7739,8328,8900,8993,9040,8467,7801,7120,6176,5325,4334,3452,3436,3405]/10;
L_h0=[1666,1643,1697,1890,2060,2515,4290,4652,4683,4891,4722,4691,4583,4699,4583,4621,4220,4328,4429,4212,3711,2237,1851,1697]/10;
Predict_wd=[4610,4423,3987,3901,3816,3481,3325,3146,3255,3325,3380,3271,3317,3200,3208,3442,3263,3107,3208,3496,3629,3652,4376,4595]/10;
pri_e=[0.40*ones(1,7),0.75*ones(1,4),1.20*ones(1,3),0.75*ones(1,4),1.20*ones(1,4),0.40*ones(1,2)];
grid_sw=[0.2*ones(1,24)];
eta_charge = 0.95;
eta_discharge = 0.96;
E_bat_initial = 800;
E_bat_max = 1800; % 最大容量
E_bat_min = 500;  % 最小容量
P_bat_charge_max = 500; % 最大总充电功率
P_bat_discharge_max = 500; % 最大总放电功率
T = 24;
%% 1. 决策变量
L_e=sdpvar(1,24); % 需求响应后电负荷
L_h=sdpvar(1,24); %微网经过需求响应后实际的热负荷
P_e_cut=sdpvar(1,24);   %微网的可削减电负荷
P_e_tran=sdpvar(1,24);  %微网的可转移电负荷
P_h_DR=sdpvar(1,24);    %微网的可削减热负荷
E_bat=sdpvar(1,24);     %微网中的储电设备的储电余量
P_batc=[112.800000000000	91	58.2000000000000	37.2000000000001	22.5000000000000	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	20.0000000000000	94	119];    %储电设备的充电功率
P_batd=[0	0	0	0	0	0	0	0	0	0	0	137.317534754402	53.4427247451345	65.6323354958296	113.075885634846	135.665528266914	0.752391102873162	0	0	0	0	0	0	0];    %储电设备的放电功率
P_e_GT=sdpvar(1,24); %燃气轮机的发电功率
P_h_GT=sdpvar(1,24); %燃气轮机的产热功率
P_h_GB=sdpvar(1,24); %余热锅炉的产热功率
P_buy=sdpvar(1,24);
P_sell=sdpvar(1,24);
Gas=sdpvar(1,24);
P_e_wd = sdpvar(1,24); % <<<<< 新增：风电实际出力决策变量

Gas_GT=sdpvar(1,24); %GT的耗气量
Gas_GB=sdpvar(1,24); %GB的耗气量

% 向SESO出租储能的变量
P_mg1_lease_charge = sdpvar(1, 24); % 出租给SESO的充电功率
P_mg1_lease_discharge = sdpvar(1, 24); % 出租给SESO的放电功率
u_total_charge = binvar(1, 24);  % 电池总充电状态 (1=充电, 0=不充电)
u_total_discharge = binvar(1, 24); % 电池总放电状态 (1=放电, 0=不放电)
lease_revenue = 0.1; % 向SESO出租储能的收益单价 (元/kWh)

%% 3. 约束条件
C = [];

for t=1:24
    C=[C,
       L_e(t)==L_e0(t)+P_e_cut(t)+P_e_tran(t),
       L_h(t)==L_h0(t)-P_h_DR(t),
       -0.15*L_e0(t)<=P_e_cut(t)<=0,
       -0.15*L_e0(t)<=P_e_tran(t)<=0.15*L_e0(t),
       0<=P_h_DR(t)<=0.2*L_h0(t),
      ];
end
C=[C,sum(P_e_tran)==0,];
for t=1:24
    C=[C,
       P_h_GT(t)==(1-0.35)/0.35*0.83*P_e_GT(t),
       0<=P_e_GT(t)<=5000,
       0<=P_h_GB(t)<=800,
       0<=P_e_wd(t)<=Predict_wd(t), % <<<<< 新增：实际出力小于等于预测值
      ];
end
for t=1:24
    C=[C,
             P_h_GT(t)+P_h_GB(t)==L_h(t),
      ];
end
for t=1:24
    C=[C,
       500>=P_buy(t)>=0,
       0<=P_sell(t)<=2000,
       0<=P_mg1_lease_charge(t)<=500,
       0<=P_mg1_lease_discharge(t)<=500,
      ];
end
for t=1:24
    C=[C,
       P_e_GT(t)==0.35*9.7*Gas_GT(t),
       P_h_GB(t)==0.9*9.7*Gas_GB(t),
       Gas(t)==Gas_GT(t)+Gas_GB(t),
      ];
end
E_co2=0.55*sum(P_e_GT)+0.65*sum(P_h_GB);
E_0=0.424*sum(P_e_GT+P_e_wd+P_h_GB); % <<<<< 修正：使用决策变量 P_e_wd
C4=0.75*(E_co2-E_0);
C3=0.01329*sum(P_e_GT);
% (省略了未改变的储能约束...)
for t = 1:24
    P_total_charge = P_batc(t) + P_mg1_lease_charge(t);
    P_total_discharge = P_batd(t) + P_mg1_lease_discharge(t);
    C = [C, u_total_charge(t) + u_total_discharge(t) <= 1];
    C = [C, 0 <= P_total_charge <= u_total_charge(t) * P_bat_charge_max];
    C = [C, 0 <= P_total_discharge <= u_total_discharge(t) * P_bat_discharge_max];
    C = [C, P_batc(t) >= 0, P_batd(t) >= 0];
    C = [C, P_mg1_lease_charge(t) >= 0, P_mg1_lease_discharge(t) >= 0];
end
P_total_charge_t1 = P_batc(1) + P_mg1_lease_charge(1);
P_total_discharge_t1 = P_batd(1) + P_mg1_lease_discharge(1);
C = [C, E_bat(1) == E_bat_initial + P_total_charge_t1 * eta_charge - P_total_discharge_t1 / eta_discharge];
for t = 2:T
    P_total_charge_t = P_batc(t) + P_mg1_lease_charge(t);
    P_total_discharge_t = P_batd(t) + P_mg1_lease_discharge(t);
    C = [C, E_bat(t) == E_bat(t-1) + P_total_charge_t * eta_charge - P_total_discharge_t / eta_discharge];
end
for t = 1:T
    C = [C, E_bat_min <= E_bat(t) <= E_bat_max];
end
C = [C, E_bat(T) == E_bat_initial];

% 3.2 电力平衡约束
for t=1:24
    total_gen = P_e_GT(t) + P_e_wd(t) + P_buy(t) + P_batd(t); % <<<<< 修正：使用决策变量 P_e_wd
    total_cons = L_e(t) + P_sell(t) + P_batc(t) ;
    C = [C, total_gen == total_cons];
end

%% 4. 目标函数 (省略了未改变的目标函数...)
Obj_original = sum(pri_e.*P_buy)+3.5*sum(Gas)+0.01*sum(abs(P_e_tran))-0.03*sum(P_e_cut)-sum(grid_sw.*P_sell)+0.01*sum(P_total_charge_t+P_total_discharge_t)+C3+C4;
Obj_lease_revenue = sum(lease_revenue * (P_mg1_lease_charge + P_mg1_lease_discharge));
admm_penalty_lease = sum(lambda_mg1_c .* (P_mg1_lease_charge - P_seso_from_mg1_c)) ...
                   + sum(lambda_mg1_d .* (P_mg1_lease_discharge - P_seso_from_mg1_d)) ...
                   + (rho/2) * (norm(P_mg1_lease_charge - P_seso_from_mg1_c)^2 + norm(P_mg1_lease_discharge - P_seso_from_mg1_d)^2);
obj_mg1 = Obj_original + Obj_lease_revenue + admm_penalty_lease ;

%% 5. 求解与输出 (省略了未改变的求解部分...)
ops = sdpsettings('solver','gurobi','verbose',0,'debug', 0);
solution = solvesdp(C, obj_mg1, ops);
if solution.problem ~= 0
    fprintf('!!! 警告: 在 Fun_MG1_Lease 中求解失败 !!!\n');
    disp(yalmiperror(solution.problem));
end
P_mg1_net_exchange = value(P_sell) - value(P_buy);
P_mg1_lease_charge = double(P_mg1_lease_charge);
P_mg1_lease_discharge = double(P_mg1_lease_discharge);
obj_mg1 = double(Obj_original);
E_bat_val = value(E_bat);
P_batc_val = value(P_batc);
P_batd_val = value(P_batd);
P_e_GT_val = value(P_e_GT);
P_e_wd_val = value(P_e_wd); % <<<<< 修正：输出正确的决策变量
P_buy_val = value(P_buy);
P_sell_val = value(P_sell);
L_e_val = value(L_e);
P_e_cut_val=value(P_e_cut);
P_e_tran_val=value(P_e_tran);
P_h_DR_val=value(P_h_DR);
  end
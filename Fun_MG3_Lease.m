% ===== 完整修正版 Fun_MG3_Lease.m =====
function [P_mg3_net_exchange, P_mg3_lease_charge, P_mg3_lease_discharge, obj_mg3, E_bat_val, P_batc_val, P_batd_val, P_e_GT_val, P_e_pv_val, P_buy_val, P_sell_val, L_e_val,P_e_cut_val,P_e_tran_val,P_h_DR_val] = Fun_MG3_Lease( P_seso_from_mg3_c, P_seso_from_mg3_d, lambda_mg3_c, lambda_mg3_d, rho)
% Fun_MG_Lease: 考虑储能租赁的微网模型
%微网参数设置
% clc; clear; close all;
L_e0=[1180,1073,1196,1165,1165,1165,1196,1503,1733,2147,2193,2407,2469,2653,2852,2208,2285,2300,3220,3220,2423,1993,1779,1518]/10;
L_h0=[1494,1448,1325,1302,1317,1494,1594,1833,2080,2211,2311,2388,2303,2526,2434,2326,2164,2126,2118,2519,1841,1625,1494,1394]/10;
Predict_pv=[0,0,0,0,0,0,663,1084,1903,2277,2386,2480,2402,2168,2012,1474,998,0,0,0,0,0,0,0]/10;
pri_e=[0.40*ones(1,7),0.75*ones(1,4),1.20*ones(1,3),0.75*ones(1,4),1.20*ones(1,4),0.40*ones(1,2)];
grid_sw=[0.2*ones(1,24)];
eta_charge = 0.95;
eta_discharge = 0.96;
E_bat_initial = 800;
E_bat_max = 1800;
E_bat_min = 500;
P_bat_charge_max = 500;
P_bat_discharge_max = 500;
T=24;
%% 1. 决策变量
L_e=sdpvar(1,24);
L_h=sdpvar(1,24);
P_e_cut=sdpvar(1,24);
P_e_tran=sdpvar(1,24);
P_h_DR=sdpvar(1,24);
E_bat=sdpvar(1,24);
P_batc=[0	0	0	0	0	0	0	53.2325301204819	124.951807228916	127.750695088044	139.240685820204	16.2508432515497	0	82.5991658943465	42.3243744207599	0	0	0	0	0	0	0	0	0];
P_batd=[5.06153846153848	0	14.9526227988879	13.9763206672846	13.1978220574607	4.01153846153849	0	0	0	0	0	0	0	0	0	0	0	50.6607970342911	128.088531047266	94.6641334569046	74.0622706209454	55.1726506024097	46.9915384615385	33.9115291936979];
P_e_GT=sdpvar(1,24);
P_h_GT=sdpvar(1,24);
P_h_GB=sdpvar(1,24);
P_buy=sdpvar(1,24);
P_sell=sdpvar(1,24);
Gas=sdpvar(1,24);
P_e_pv = sdpvar(1,24); % <<<<< 修正/新增：光伏实际出力决策变量
Gas_GT=sdpvar(1,24);
Gas_GB=sdpvar(1,24);
P_mg3_lease_charge = sdpvar(1, 24);
P_mg3_lease_discharge = sdpvar(1, 24);
u_total_charge = binvar(1, 24);
u_total_discharge = binvar(1, 24);
lease_revenue = 0.1;

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
       0<=P_e_pv(t)<=Predict_pv(t), % <<<<< 新增：实际出力小于等于预测值
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
       0<=P_sell(t)<=20,
       0<=P_mg3_lease_charge(t)<=500,
       0<=P_mg3_lease_discharge(t)<=500,
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
E_0=0.424*sum(P_e_GT+P_e_pv+P_h_GB); % <<<<< 修正：使用决策变量 P_e_pv
C4=0.75*(E_co2-E_0);
C3=0.01329*sum(P_e_GT);

for t = 1:24
    P_total_charge = P_batc(t) + P_mg3_lease_charge(t);
    P_total_discharge = P_batd(t) + P_mg3_lease_discharge(t);
    C = [C, u_total_charge(t) + u_total_discharge(t) <= 1];
    C = [C, 0 <= P_total_charge <= u_total_charge(t) * P_bat_charge_max];
    C = [C, 0 <= P_total_discharge <= u_total_discharge(t) * P_bat_discharge_max];
    C = [C, P_batc(t) >= 0, P_batd(t) >= 0];
    C = [C, P_mg3_lease_charge(t) >= 0, P_mg3_lease_discharge(t) >= 0];
end
P_total_charge_t1 = P_batc(1) + P_mg3_lease_charge(1);
P_total_discharge_t1 = P_batd(1) + P_mg3_lease_discharge(1);
C = [C, E_bat(1) == E_bat_initial + P_total_charge_t1 * eta_charge - P_total_discharge_t1 / eta_discharge];
for t = 2:T
    P_total_charge_t = P_batc(t) + P_mg3_lease_charge(t);
    P_total_discharge_t = P_batd(t) + P_mg3_lease_discharge(t);
    C = [C, E_bat(t) == E_bat(t-1) + P_total_charge_t * eta_charge - P_total_discharge_t / eta_discharge];
end
for t = 1:T
    C = [C, E_bat_min <= E_bat(t) <= E_bat_max];
end
C = [C, E_bat(T) == E_bat_initial];

% 3.2 电力平衡约束
for t=1:24
    total_gen = P_e_GT(t) + P_e_pv(t) + P_buy(t) + P_batd(t); % <<<<< 修正：使用决策变量 P_e_pv
    total_cons = L_e(t) + P_sell(t) + P_batc(t) ;
    C = [C, total_gen == total_cons];
end

%% 4. 目标函数
Obj_original = sum(pri_e.*P_buy)+3.5*sum(Gas)+0.01*sum(abs(P_e_tran))-0.03*sum(P_e_cut)-sum(grid_sw.*P_sell)+0.01*sum(P_total_charge_t+P_total_discharge_t)+C3+C4;
Obj_lease_revenue = sum(lease_revenue * (P_mg3_lease_charge + P_mg3_lease_discharge));
admm_penalty_lease = sum(lambda_mg3_c .* (P_mg3_lease_charge - P_seso_from_mg3_c)) ...
                   + sum(lambda_mg3_d .* (P_mg3_lease_discharge - P_seso_from_mg3_d)) ...
                   + (rho/2) * (norm(P_mg3_lease_charge - P_seso_from_mg3_c)^2 + norm(P_mg3_lease_discharge - P_seso_from_mg3_d)^2);
obj_mg3 = Obj_original + Obj_lease_revenue + admm_penalty_lease ;

%% 5. 求解与输出
ops = sdpsettings('solver','gurobi','verbose',0,'debug', 0);
solution = solvesdp(C, obj_mg3, ops);
if solution.problem ~= 0
    fprintf('!!! 警告: 在 Fun_MG3_Lease 中求解失败 !!!\n');
    disp(yalmiperror(solution.problem));
end
% <<<<< 修正：此处函数输出变量名有误，应为 P_mg3_net_exchange
P_mg3_net_exchange = value(P_sell) - value(P_buy);
P_mg3_lease_charge = double(P_mg3_lease_charge);
P_mg3_lease_discharge = double(P_mg3_lease_discharge);
obj_mg3 = double(Obj_original);
E_bat_val = value(E_bat);
P_batc_val = value(P_batc);
P_batd_val = value(P_batd);
P_e_GT_val = value(P_e_GT);
P_e_pv_val = value(P_e_pv); % <<<<< 修正：输出正确的决策变量
P_buy_val = value(P_buy);
P_sell_val = value(P_sell);
L_e_val = value(L_e);
P_e_cut_val=value(P_e_cut);
P_e_tran_val=value(P_e_tran);
P_h_DR_val=value(P_h_DR);
end
function [P_mg2_net_exchange, P_mg2_lease_charge, P_mg2_lease_discharge, obj_mg2, E_bat_val, P_batc_val, P_batd_val, P_e_GT_val, P_e_pv_val, P_buy_val, P_sell_val, L_e_val,P_e_cut_val,P_e_tran_val,P_h_DR_val] = Fun_MG2_Lease( P_seso_from_mg2_c, P_seso_from_mg2_d, lambda_mg2_c, lambda_mg2_d, rho)
% Fun_MG_Lease: 考虑储能租赁的微网模型
%微网参数设置
%clc; clear; close all;
L_e0=[1774,1450,1296,1219,1095,1265,1481,1944,2484,2083,1651,1188,1080,1126,1033,1033,941,1450,2283,3148,3904,3719,2746,2453]/10;
L_h0=[1610,1594,1594,1610,1633,1633,1286,1201,1117,1109,1648,1656,1664,1140,1124,1109,1286,1309,1302,1325,1479,1502,1340,1332]/10;
Predict_pv=[0,0,0,0,0,0,967,1287,1583,1833,1918,1942,2004,1957,1669,1076,655,0,0,0,0,0,0,0]/10;
pri_e=[0.40*ones(1,7),0.75*ones(1,4),1.20*ones(1,3),0.75*ones(1,4),1.20*ones(1,4),0.40*ones(1,2)];
grid_sw=[0.2*ones(1,24)];
eta_charge = 0.95;
eta_discharge = 0.96;
E_bat_initial = 800;
E_bat_max = 1800; % 最大容量
E_bat_min = 500;  % 最小容量
P_bat_charge_max = 500; % 最大总充电功率
P_bat_discharge_max = 500; % 最大总放电功率
T=24;

%% 1. 决策变量
L_e=sdpvar(1,24);
L_h=sdpvar(1,24);
P_e_cut=sdpvar(1,24);
P_e_tran=sdpvar(1,24);
P_h_DR=sdpvar(1,24);
E_bat=sdpvar(1,24);
P_batc=[0	0	0	0	0	0	59.7732808155700	0	0	95.0469972196478	161.761047265987	111.040000000000	211.161445783133	116.880000000000	152.925495829472	71.6334456042798	66.3732808155699	0	0	0	0	0	0	0];
P_batd=[40.6211492122336	18.7715477293791	7.99154772937906	1.77114921223361	0	3.79745134383693	0	0	0	0	0	0	0	0	0	0	0	33.5630213160334	92.2363206672846	151.592622798888	196.520037071362	182.376339202966	122.674124189064	102.579323447637];
P_e_GT=sdpvar(1,24);
P_h_GT=sdpvar(1,24);
P_h_GB=sdpvar(1,24);
P_buy=sdpvar(1,24);
P_sell=sdpvar(1,24);
Gas=sdpvar(1,24);
Gas_GT=sdpvar(1,24);
Gas_GB=sdpvar(1,24);
P_e_pv = sdpvar(1,24); 


P_mg2_lease_charge = sdpvar(1, 24);
P_mg2_lease_discharge = sdpvar(1, 24);
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
for t=1:24
    C=[C,
       P_h_GT(t)==(1-0.35)/0.35*0.83*P_e_GT(t),
       0<=P_e_GT(t)<=5000,
       0<=P_h_GB(t)<=800,
       0<=P_e_pv(t)<=Predict_pv(t), 
      ];
end
for t=1:24
    C=[C,
             P_h_GT(t)+P_h_GB(t)==L_h(t),
      ];
end
for t=1:24
    C=[C,
       2000>=P_buy(t)>=0,
       0<=P_sell(t)<=2000,
       0<=P_mg2_lease_charge(t)<=500,
       0<=P_mg2_lease_discharge(t)<=500,
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
E_0=0.424*sum(P_e_GT + P_e_pv + P_h_GB); 
C4=0.75*(E_co2-E_0);
C3=0.01329*sum(P_e_GT);

% --- 储能约束 (无修改) ---
for t = 1:24
    P_total_charge = P_batc(t) + P_mg2_lease_charge(t);
    P_total_discharge = P_batd(t) + P_mg2_lease_discharge(t);
    C = [C, u_total_charge(t) + u_total_discharge(t) <= 1];
    C = [C, 0 <= P_total_charge <= u_total_charge(t) * P_bat_charge_max];
    C = [C, 0 <= P_total_discharge <= u_total_discharge(t) * P_bat_discharge_max];
    C = [C, P_batc(t) >= 0, P_batd(t) >= 0];
    C = [C, P_mg2_lease_charge(t) >= 0, P_mg2_lease_discharge(t) >= 0];
end
P_total_charge_t1 = P_batc(1) + P_mg2_lease_charge(1);
P_total_discharge_t1 = P_batd(1) + P_mg2_lease_discharge(1);
C = [C, E_bat(1) == E_bat_initial + P_total_charge_t1 * eta_charge - P_total_discharge_t1 / eta_discharge];
for t = 2:T
    P_total_charge_t = P_batc(t) + P_mg2_lease_charge(t);
    P_total_discharge_t = P_batd(t) + P_mg2_lease_discharge(t);
    C = [C, E_bat(t) == E_bat(t-1) + P_total_charge_t * eta_charge - P_total_discharge_t / eta_discharge];
end
for t = 1:T
    C = [C, E_bat_min <= E_bat(t) <= E_bat_max];
end
C = [C, E_bat(T) == E_bat_initial];


% --- 3.2 电力平衡约束 ---
for t=1:24
    total_gen = P_e_GT(t) + P_e_pv(t) + P_buy(t) + P_batd(t); % ★★★★★ 关键修正 4: 在平衡中使用决策变量
    total_cons = L_e(t) + P_sell(t) + P_batc(t) ;
    C = [C, total_gen == total_cons];
end


%% 4. 目标函数
Obj_original = sum(pri_e.*P_buy)+3.5*sum(Gas)+0.01*sum(abs(P_e_tran))-0.03*sum(P_e_cut)-sum(grid_sw.*P_sell)+0.01*sum(P_total_charge_t+P_total_discharge_t)+C3+C4; 
Obj_lease_revenue = sum(lease_revenue * (P_mg2_lease_charge + P_mg2_lease_discharge));
admm_penalty_lease = sum(lambda_mg2_c .* (P_mg2_lease_charge - P_seso_from_mg2_c)) ...
                   + sum(lambda_mg2_d .* (P_mg2_lease_discharge - P_seso_from_mg2_d)) ...
                   + (rho/2) * (norm(P_mg2_lease_charge - P_seso_from_mg2_c)^2 + norm(P_mg2_lease_discharge - P_seso_from_mg2_d)^2);
obj_mg2 = Obj_original + Obj_lease_revenue + admm_penalty_lease ;

%% 5. 求解与输出
ops = sdpsettings('solver','gurobi','verbose',0,'debug', 0);
solution = solvesdp(C, obj_mg2, ops);
if solution.problem ~= 0
    fprintf('!!! 警告: 在 Fun_MG2_Lease 中求解失败 !!!\n');
    disp(yalmiperror(solution.problem));
end
P_mg2_net_exchange = value(P_sell) - value(P_buy);
P_mg2_lease_charge = double(P_mg2_lease_charge);
P_mg2_lease_discharge = double(P_mg2_lease_discharge);
obj_mg2 = double(Obj_original);
E_bat_val = value(E_bat);
P_batc_val = value(P_batc);
P_batd_val = value(P_batd);
P_e_GT_val = value(P_e_GT);
P_e_pv_val = value(P_e_pv); 
P_buy_val = value(P_buy);
P_sell_val = value(P_sell);
L_e_val = value(L_e);
P_e_cut_val=value(P_e_cut);
P_e_tran_val=value(P_e_tran);
P_h_DR_val=value(P_h_DR);
 end
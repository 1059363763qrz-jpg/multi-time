    function [P_dso_charge, P_dso_discharge, obj_dso, P_grid_exchange, PLOAD_val, P_loss_val, P_pv_val, P_wt_val, PLOAD_original_val, Lshift_effect_val, SIL_effect_val] = Fun_DSO( P_seso_to_dso_c, P_seso_to_dso_d, lambda_dso_c, lambda_dso_d, P_mg1_net_exchange, P_mg2_net_exchange, P_mg3_net_exchange, rho)
   % clc; clear; close all;
%% 1.设参
mpc = IEEE33BW;
pload = mpc.Pload*100;%节点有功负荷,原数据100MVA，改为1MVA
PPload=[1.77184023	1.760691425	1.985343246	2.164759614	2.525216573	2.836282757	2.717289933	2.912703959	2.938170633	3.200578131	3.93626314	3.680872774	3.53168293	3.215685929	2.903666471	2.626416093	2.33835944	2.582209352	3.063339069	3.465731913	3.01938202	2.591419695	2.254491338	1.784180138]*5;
qload = mpc.Qload*100;%节点无功负荷
branch = mpc.branch;
branch(:,3) = branch(:,3)*1/(10^2);%求阻抗标幺值
R = real(branch(:,3));
X = imag(branch(:,3));
T = 24;%时段数为24小时
nb = 33;%节点数
nl = 32;%支路数
nwt = 2;%风机
npv = 1;%光伏
S_pv=1.5;
RR=xlsread('光照强度数据.xlsx','B2:B25');
pp_pv=xlsread('光照强度数据.xlsx','D3:AA3');
C_e=xlsread('电价.xlsx','A1:X1');
C_ae=xlsread('激励电价.xlsx','A1:X1');
C_de=xlsread('需求响应激励价格.xlsx','A1:X1');
R_STD=1000;%标准太阳辐射
R_C=150;%辐射点
PP_pv=zeros(npv,T);
% %%%%光伏出力模型
for i=1:1:npv 
    for t=1:1:T
        if(0<=RR(t))&(RR(t)<=R_C)
            PP_pv(i,t)=S_pv*((RR(t))^2)/(R_STD*R_C);
        else if(R_C<=RR(t))&(RR(t)<=R_STD)
                PP_pv(i,t)=S_pv*RR(t)/R_STD;
            else if(R_STD<=RR(t))
                    PP_pv(i,t)=S_pv;
                end
            end
        end
    end
end
v=xlsread('风速数据.xlsx','B2:B25');
pp_wt=xlsread('风速数据.xlsx','G2:AD2');
v_in=3;%切入风速
v_r=11.3;%额定风速
v_out=25;%切出风速
S_wt=[1;1.5];
PP_wt=zeros(nwt,T);
%%%风电出力模型
for i=1:1:nwt 
    for t=1:1:T
      if (0<=v(t)) &(v(t)<=v_in)|(v_out<=v(t))
          PP_wt(i,t)=0;
      else if(v_in<=v(t))&(v(t)<v_r)
              
              PP_wt(i,t)=[(v(t)-v_in)/(v_r-v_in)]*S_wt(i);
          else if(v_r<=v(t))&(v(t)<v_out)
                  PP_wt(i,t)=S_wt(i);
              end
          end
      end
    end
end
ness = 1;%ESS数
v_in=3;%切入风速
v_r=11.3;%额定风速
v_out=25;%切出风速
upstream = zeros(nb,nl);
dnstream = zeros(nb,nl);
for t = 1:nl
    upstream(t,t)=1;
end
for t = [1:16,18:20,22:23,25:31]
    dnstream(t,t+1)=1;
end
dnstream(1,18) = 1;
dnstream(2,22) = 1;
dnstream(5,25) = 1;
dnstream(33,1) = 1;
Vmax = [1.1*1.1*ones(nb-1,T)
        1.1*1.1*ones(1,T)];
Vmin = [0.9*0.9*ones(nb-1,T)
        0.9*0.9*ones(1,T)];
Pgmax = [zeros(nb-1,T)
        50*ones(1,T)];
Pgmin = [zeros(nb-1,T)
         -20*ones(1,T)]; 
Qgmax = [zeros(nb-1,T)
         30*ones(1,T)];
Qgmin = [zeros(nb-1,T)
         -10*ones(1,T)];
%% 2.设变量
V = sdpvar(nb,T);%电压的平方
I = sdpvar(nl,T);%电流的平方
P = sdpvar(nl,T);%线路有功
Q = sdpvar(nl,T);%线路无功
Pg = sdpvar(nb,T);%发电机有功
Qg = sdpvar(nb,T);%发电机无功
p_wt = sdpvar(nwt,T);%风机有功
p_pv = sdpvar(npv,T);%光伏有功
s_IL=binvar(1,T);
Temp_shift=binvar(3,T);
Lshift=sdpvar(3,T);
Lshift_old=xlsread('可转移负荷.xlsx','A1:X3');
PLOAD=sdpvar(nb,T);
Pload=sum(PLOAD);
P_dso_discharge = sdpvar(ness,T);%ESS放电功率
P_dso_charge = sdpvar(ness,T);%ESS充电功率
u_dch = binvar(ness,T);%ESS放电状态
u_ch = binvar(ness,T);%ESS充电状态
E_ess = sdpvar(ness,T);%ESS的电量
S_IL1=sdpvar(1,T,'full');
S_IL2=sdpvar(1,T,'full');
%% 3.设约束
C = [];
%% 需求响应约束
%可平移负荷1约束,9-14为可平移负荷区间
 C = [C,sum(Temp_shift(1,1:8))==0,sum(Temp_shift(1,9:14))==1,sum(Temp_shift(1,15:24))==0];
 %可平移负荷2约束,19-22为可平移负荷区间
 C = [C,sum(Temp_shift(2,1:18))==0,sum(Temp_shift(2,19:22))==1,sum(Temp_shift(2,23:24))==0];
  %可平移负荷3约束,8-21为可平移负荷区间
 C = [C,sum(Temp_shift(3,1:7))==0,sum(Temp_shift(3,8:21))==1,sum(Temp_shift(3,22:24))==0];
 %可平移负荷1&可平移负荷2&可平移负荷3平移
 for k = 1:24
  C = [C, implies(Temp_shift(1,k),Lshift(1,:)==circshift(Lshift_old(1,:),[0,k-17]))];
  C = [C, implies(Temp_shift(2,k),Lshift(2,:)==circshift(Lshift_old(2,:),[0,k-3]))];
  C = [C, implies(Temp_shift(3,k),Lshift(3,:)==circshift(Lshift_old(3,:),[0,k-23]))];
 end
 TL1=[zeros(3,T);Lshift_old(1,:)-Lshift(1,:);zeros(29,T)];    
 TL2=[zeros(14,T);Lshift_old(2,:)-Lshift(2,:);zeros(18,T)];
 TL3=[zeros(29,T);Lshift_old(3,:)-Lshift(3,:);zeros(3,T)];
%%%%可削减负荷%%%%%
C=[C,0<=S_IL1<=0.8*pload(10,:)];
C=[C,0<=S_IL2<=0.8*pload(26,:)];
IL1=[zeros(9,T);S_IL1(1,:);zeros(23,T)];
IL2=[zeros(25,T);S_IL2(1,:);zeros(7,T)];
C=[C,pload+TL1+TL2+TL3-IL1-IL2==PLOAD];

%% 储能装置（ESS）约束
%充放电状态约束
C = [C, u_dch(1,:) + u_ch(1,:) <= 1];%表示充电，放电，不充不放三种状态

%充放电时刻约束
% C = [C, [sum(u_dch(1,8:15))==0,sum(u_dch(1,18:21))==0]];
% C = [C, [sum(u_ch(1,1:7))==0,sum(u_ch(1,14:18))==0,sum(u_ch(1,22:24))==0]];

%功率约束
C = [C, 0 <= P_dso_discharge <= u_dch*1];
C = [C, 0 <= P_dso_charge<= u_ch*1];
%容量约束
% for t = 1:23
%         C = [C, E_ess(1,t+1) == E_ess(1,t) + 0.9*p_ch(1,t) - 0.9*p_dch(1,t)]; 
% 
% end
% C = [C, E_ess(1,1) == E_ess(1,24) ]; 
% C = [C,  0.1<=E_ess(1,:)/1.03<=0.9];

%投入节点选择
P_dch1 = [zeros(16,T);P_dso_discharge(1,:);zeros(16,T)];
P_ch1 = [zeros(16,T);P_dso_charge(1,:);zeros(16,T)];
%与微网交互量
P_mg1_net_exchange = P_mg1_net_exchange/1000;
P_mg2_net_exchange = P_mg2_net_exchange/1000;
P_mg3_net_exchange = P_mg3_net_exchange/1000;
P_exchange = [zeros(2,T);P_mg1_net_exchange(1,:);zeros(10,T);P_mg2_net_exchange(1,:);zeros(10,T);P_mg3_net_exchange(1,:);zeros(8,T)];
%与seso交互，转化成mw
P_seso_to_dso_c=P_seso_to_dso_c/1000;
P_seso_to_dso_d=P_seso_to_dso_d/1000;
%% 风机+光伏约束
C = [C, 0 <= p_wt,p_wt <= 1.2*PP_wt];
P_wt = [zeros(15,T);p_wt(1,:);zeros(10,T);p_wt(2,:);zeros(6,T)];
C = [C, 0 <= p_pv,p_pv <=1.2*PP_pv];
p_pv =1.2*PP_pv;
P_pv = [zeros(15,T);p_pv;zeros(17,T)];
%% 潮流约束
Pin = -upstream*P + upstream*(I.*(R*ones(1,T))) + dnstream*P;%节点注入有功
Qin = -upstream*Q + upstream*(I.*(X*ones(1,T))) + dnstream*Q;%节点注入无功
   C = [C, Pin+pload+TL1+TL2+TL3-IL1-IL2- P_exchange-P_wt- P_pv- Pg-P_dch1+ P_ch1==0];
  % C = [C, Pin+pload+TL1+TL2+TL3-IL1-IL2-P_wt- P_pv- Pg-P_dch1+ P_ch1==0];
C = [C, Qin + qload - Qg == 0];
%欧姆定律约束
C = [C, V(branch(:,2),:) == V(branch(:,1),:) - 2*(R*ones(1,24)).*P - 2*(X*ones(1,24)).*Q + ((R.^2 + X.^2)*ones(1,24)).*I];
%二阶锥约束
C = [C, V(branch(:,1),:).*I >= P.^2 + Q.^2];
%% 通用约束
%节点电压约束
C = [C, Vmin <= V,V <= Vmax];
%发电机功率约束
C = [C, Pgmin <= Pg,Pg <= Pgmax,Qgmin <= Qg,Qg <= Qgmax];
%支路电流约束
C = [C, 0 <= I,I <= 1000];
%% 4.设目标函数
Pload2=PPload+sum(P_dso_discharge-P_dso_charge);
Pload3=Pload+sum(P_dso_discharge-P_dso_charge);
%%%%%%%%%
for t=1:1:23
    f1_0=0.25*[(max(PPload)-min(PPload))/mean(PPload(:))]+0.75*[(max(abs(PPload(t)-Pload(t+1))))/mean(PPload(:))];
     f1_1=0.25*[(max(Pload)-min(Pload))/mean(Pload(:))]+0.75*[(max(abs(Pload(t)-Pload(t+1))))/mean(Pload(:))];
     f1_2=0.25*[(max(Pload2)-min(Pload2))/mean(Pload2)]+0.75*[(max(abs(Pload2(t)-Pload2(t+1))))/mean(Pload2(:))];
     f1_3=0.25*[(max(Pload3)-min(Pload3))/mean(Pload3)]+0.75*[(max(abs(Pload3(t)-Pload3(t+1))))/mean(Pload3(:))];
end
%%%%%%%电压差最小%%%%%%%%%
f1=0;
for t=1:1:T
    for n=1:1:nb-1
        f1=f1+abs(V(n,t)-V(n+1,t));
    end
end
%%%%%%%%%
% B1=sum(sum(p_dso_discharge-p_dso_charge).*C_ae);
% B2=sum(sum(p_dso_discharge))*200;
% R_ess=B1+B2;   %储能补贴和收益%%%%%
C_LA=sum([0.062,0.043,0.048]*Temp_shift.*C_de)+sum((S_IL1+S_IL2).*C_de);%%%负荷调用成本%%%%
C_buy=sum(sum(Pg).*C_e);
C_loss=sum(sum(I.*(R*ones(1,T))))*400;  %%%网损成本%%%% 
% C_ess=250*sum(sum(p_dso_discharge)+sum(p_dso_charge));  %%%储能调用成本%%%%
 C_lease=250*sum(P_dso_discharge)-550*sum(P_dso_charge);%优先用储能
f2=C_LA+C_loss+C_buy;

admm_penalty = sum(lambda_dso_c .* (P_dso_charge - P_seso_to_dso_c)) ...
             + sum(lambda_dso_d .* (P_dso_discharge - P_seso_to_dso_d)) ...
             + (rho*1000/2) * (norm(P_dso_charge - P_seso_to_dso_c)^2 + norm(P_dso_discharge - P_seso_to_dso_d)^2);
obj_dso = f2+C_lease+admm_penalty;
%%%%储能寿命损耗成本%%%%
% n=sum(sum(u_dch + u_ch));
% L=0;
% for D=0:0.1:1 
%     L=L+n/(71470*D^4-170100*D^3+146400*D^2-56500*D+12230);
% end
% C_ESS1=819000*2+2953000*0.375;
% C_ESS2=0;
% for t=1:1:10
%     C_ESS2=C_ESS2+69000*2*((1+0.015)/(1+0.09))^t;
% end
% C_ESS=C_ESS1+C_ESS2;
% C_day=C_ESS*L;       %%%%%储能寿命损耗成本%%%%%%
%%%%%%%%%%%%%%%
%%%%%%配电网综合收益%%%%%%
peak_0=PPload(9)+PPload(10)+PPload(11)+PPload(12)+PPload(13)+PPload(14)+PPload(19)+PPload(20);
peak_3=Pload3(9)+Pload3(10)+Pload3(11)+Pload3(12)+Pload3(13)+Pload3(14)+Pload3(19)+Pload3(20);
lamda3=(peak_0-peak_3)/peak_0;  %%%%场景三的削峰率
C_gridup=1;
delta_n3=((log10(1+lamda3))/log10(1+0.015));   %%%%%%延缓改造年限%%%%
F1=C_gridup*[1-((1+0.015)/(1+0.09))^delta_n3];  %%%%%%减少配电网升级改造费用%%%%
T=9;   %%%%充放电周期%%%%
e_s=1000;   %%%%发电厂度电收益%%%%%
F2=0.5*0.375*T*e_s;   %%%%%减少备用成本%%%%%%
%toc%建模时间
%% 5.设求解器
warning('off', 'YALMIP:bigm'); 
ops = sdpsettings('verbose', 0, 'solver', 'gurobi');
sol = optimize(C,obj_dso,ops);
if sol.problem == 0
    disp('DSO subproblem solved successfully');
P=value(P)*1000;
 P_dso_charge = value(P_dso_charge) * 1000;
    P_dso_discharge = value(P_dso_discharge) * 1000;
    obj_dso = value(f2);
    
    % P_buy_val 是与主网的净交换功率 (购电为正, 售电为负)
    P_grid_exchange = value(Pg(33,:)) * 1000;
    
    PLOAD_val = value(PLOAD) * 1000;
    P_pv_val = value(p_pv) * 1000;
    P_wt_val = value(p_wt) * 1000;
    
    % P_exchange_val 是与微网的净交换功率 (DSO购电为正, 售电为负)
    P_exchange_val = value(P_exchange) * 1000;

    % 步骤 2: 使用功率平衡法计算网损
    % 定义总电源 (Source)，包括所有注入DSO的功率
    total_source = P_grid_exchange + sum(P_pv_val, 1) + sum(P_wt_val, 1) + P_dso_discharge ...
                   + sum(P_exchange_val, 1); 
                   
    % 定义总负荷 (Load)，包括所有从DSO流出的功率
    total_load = sum(PLOAD_val, 1) + P_dso_charge;
                 
    % 计算网损 (损耗 = 注入 - 流出)
    P_loss_val = total_source - total_load;
    
    % --- 为保证函数输出完整，补充定义其他返回值 ---
    PLOAD_original_val = pload * 1000;
    Lshift_effect_val = value(sum(TL1 + TL2 + TL3, 1)) * 1000;
    SIL_effect_val = value(sum(IL1 + IL2, 1)) * 1000;
    % ======================== END OF DEFINITIVE VERSION ========================

else
    disp('!!! DSO subproblem FAILED to solve !!!');
    yalmiperror(sol.problem);
    
    % 当求解失败时，将所有输出变量设为 NaN
    P_dso_charge = nan(ness,T); P_dso_discharge = nan(ness,T); obj_dso = NaN;
    P_grid_exchange = nan(1,T); PLOAD_val = nan(nb,T); P_loss_val = nan(1,T);
    P_pv_val = nan(npv,T); P_wt_val = nan(nwt,T); PLOAD_original_val = pload*1000;
    Lshift_effect_val = nan(1,T); SIL_effect_val = nan(1,T);

     end
%toc%求解时间
clear branch C dnstream upstream i kk mpc nb nl ncb ness noltc npv nwt nsvc...
      QCB_step ops Pgmax Pgmin pload qload t T theta_DE theta_IN u_ch u_dch...
      Vmax Vmin R X P_ch P_dch P_pv P_wt Q_SVC Qgmax Qgmin Pin Qin theta_CB...
      Q_CB k r rjs theta1_DE theta1_IN theta_OLTC
%% 6.分析错误标志
if sol.problem == 0
    disp('succcessful solved');
else
    disp('error');
    yalmiperror(sol.problem)
end
% %% 绘图
% 
% 
% figure(3)
% t=1:1:24;
% plot(t,sum(p_dch-p_ch),'-^k','linewidth',1.8);
% hold on
% plot(t,Pg(33,:),'-dm','linewidth',1.8);
% xlabel('时刻/h');
% ylabel('有功功率/MW');
% legend('储能充放电功率','主网出力');
% 
% figure(4)
% t=1:1:24;
% plot(t,PPload,'-','linewidth',1.5);
% hold on
% plot(t,PPload+sum(Lshift_old-Lshift)-S_IL1-S_IL2,'-.','linewidth',1.5);
% hold on
% plot(t,PPload+sum(p_dch-p_ch),'-d','linewidth',1.5);
% hold on
% plot(t,PPload+sum(p_dch-p_ch)+sum(Lshift_old-Lshift)-S_IL1-S_IL2,'-*','linewidth',1.5);
% xlabel('时刻/h');
% ylabel('负荷有功功率/MW')
% legend('实际值','场景一','场景二','场景三');
% 
% figure(5)
% t=1:1:24;
% plot(t,p_dch(1,:)-p_ch(1,:),'-.','linewidth',1.8);
% hold on
% 
% xlabel('时刻/h');
% ylabel('有功功率/MW');
% legend('ESS1');
% 
% figure(6)
% t=1:1:24;
% plot(t,E_ess(1,:)/1.03,'-.','linewidth',1.8);
% hold on
% 
% xlabel('时刻/h');
% ylabel('SOC');
% legend('ESS1');
% 
% 
% 
% figure(10)
% t=1:1:24;
% S_DR=zeros(5,24);
% S_DR(1:3,:)=Lshift_old-Lshift;
% S_DR(4,:)=-S_IL1;
% S_DR(5,:)=-S_IL2;
% bar(t,(S_DR)','stacked');
% xlabel('时刻/h');
% ylabel('负荷功率/MW');
% legend('可转移负荷1','可转移负荷2','可转移负荷3','可削减负荷1','可削减负荷2');
     % end
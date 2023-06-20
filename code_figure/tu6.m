%---------------------------------基本信息---------------------------------%
% 程序名称：bpsk_ber_v6
% 作者信息：honney_li
% 实现功能：计算加入大气、混沌、相位噪声情况下的BPSK误码率
% 编译环境：matlab2018b
% 编写日期：2022/10/12
%-------------------------------------------------------------------------%

%-------------------------------版本说明-------------------------------%

%{

1.在V2版本的基础上加入混沌噪声的BER   
2.分析了4种激光器的线宽，分别为500 800 1200 2000KHz
3.在V3版本的基础上加入 强湍流 弱湍流 还有夜视情况下的大气湍流 PDF 两重积分 单独一重pdf不等于0
4.在V4版本的基础上加入重新整理的公式部分
5.具体版本说明的word 下载网址：
https://artic-1302263000.cos.ap-nanjing.myqcloud.com/BER%E6%8E%A8%E5%AF%BC%E6%94%B93.docx
6.重新梳理相关公式，根据公式修改程序

%}

%-------------------------------版本问题-------------------------------%
%{

1.beta的取值问题，曾学长公式是定值，但是他的代码是积分形式 曾学长小论文公式11   
2.SNR的影响问题，噪声的函数和光概率密度函数不匹配
3.算失配噪声的时候，没有用到掩膜比，这个物理量的作用需要问下
4.要分析更多激光器的线宽，将线宽也分部积分，不能只分析4种

%}

%-------------------------------程序主体-------------------------------%
clear all;
close all;
clc;
format long;%有效数字16位

%-------------------------------系统参数-------------------------------%
%一些物理常量
hp = 6.6260693 * 10^(-34);          %普朗克常量 
e = 1.60217733 * 10^(-19);          %电子电量
kc = 1.3806505 * 10^(-23);          %玻尔兹曼常数
c = 3 * 10^8;                       %光速

%-------------------------------链路相关参数-------------------------------%
zenith_angle = 0;                   %天顶角，单位：°
zeta = deg2rad(zenith_angle);       %ζ天顶角，单位：rad
H = 38.*10.^6;                      %是卫星距离地面的高度，单位:m
h0 = 100;                           %在下行中是接收站高度，在上行中是发射端高度，单位：m
L=(H-h0).*sec(zeta);                %链路距离，单位：m
alpha = 1;                          %α链路的能量损耗系数，为1代表没有损耗

%-------------------------------发射端参数-------------------------------%
P1 = 8*1e-3;                       %发送端激光器恒定功率
lambda = 1550.*10^(-9);             %λ波长，单位：m
k = 2 .* pi ./ lambda;              %波数，单位：rad/m
v = c ./ lambda;                    %ν= c/λ，光频率，单位：Hz
%divergence_angle = 1;%束散角（全角），单位：°
%theta = deg2rad(divergence_angle);%θ束散角（全角），单位：rad
theta = 30  .* 10^(-6);             %θ束散角（全角），单位：rad
W0 = 0.1;                           %发射口径,出射光斑大小，单位：m，保证光斑半径大于束腰半径
GEDFA1 = 2000;                      %发射端的EDFA增益
maskeff = 1.2;                      %掩膜比 信号功率比上载波功率

%-------------------------------调制端参数-------------------------------%
Ps = maskeff * P1 ;                       %调制端 信号激光器恒定功率 

%-------------------------------接收端参数-------------------------------%
P2 = P1/2;                        %接收端激光器恒定功率,发射功率的一半，为了混沌相消
W = W0 + theta .* L ./ 2;          %接受面上的光斑半径 W = θL/2，为近似公式
Dr = 0.40;                         %接收口径,单位：m
Waist=2 .* lambda ./ (pi.*theta);  %束腰半径 
f = pi * (Waist.^2) / lambda;      %共焦参数
Rr = L + f .^ 2 ./ L;              %接收面上的波前曲率半径
%GEDFA2 = 400;                     %接收端的EDFA增益

                                     
%-------------------------------BPSK解调模块相关参数-------------------------------%

%-------------------------------失配噪声参数-------------------------------%

%beta=10;                           %chaos feedback gain 环路增益

%-------------------------------探测器参数-------------------------------%
%%平衡探测器参数
yita=0.8;                        %η量子效率
R_d=yita*e/(hp*v);                %计算PD响应度R_d

% 探测器噪声参数设置
G=1;                             %探测器增益
F_A=G^0.5;
I_d=10^-9;                       %暗电流
delt_F=2*10^9;                 %噪声等效带宽
Ts=1/delt_F;                     %噪声等效时宽
T=300;                           %探测器温度
R_L=50;                          %探测器负载电阻
k_B=1.380649*10^-23;             %波尔兹曼常数
F_n=1;                           %探测器噪声因子


%-------------------------------  大气链路光强计算  -------------------------------%
%r=0接收光脉冲平均功率大小 经过大气链路损耗，接收端中心光功率
Pt = (P1 + Ps);                                    %混沌掩盖比的光功率
r = 0;                                               %不考虑对准误差
I_0L = alpha .* Pt .*GEDFA1.* Dr.^2 ./ (2 .* W.^2);  %大气链路的损耗,经过EDFA2之前
 
alphabeta = Dr.^2 ./ (2 .* W.^2);

%-------------------------------  大气光强闪烁参数计算   -------------------------------%


sigma_I_0L_dark =0.05;
sigma_I_0L_ruo=0.10; %sigma_I_0L大气光强闪烁方差 弱湍流

sigma_I_0L_qiang1=0.15; %sigma_I_0L大气光强闪烁方差 强湍流
sigma_I_0L_qiang2=0.20;
%大气影响下的光强概率密度函数
num = 1000+1;

I = linspace(10^-15,3*10^-5,num);               %接收光脉冲功率大小,为变量  为了最后加权计算BER
dI =(max(I)-min(I))./num;                      %步长

%-------------------------------  下行链路概率密度函数计算   -------------------------------%
%强湍流情况下

Pr_I_qiang_test = 1./(sqrt(2*pi.*sigma_I_0L_qiang1))./(I).*exp(-(log(I./I_0L) + ...
            sigma_I_0L_qiang1./2).^2./(2.*sigma_I_0L_qiang1));
       
Pr_I_qiang1 = (2 .* pi .* sigma_I_0L_qiang1 .* I .* I).^(-0.5) .* ...
                 exp(-(log(I ./ I_0L) + sigma_I_0L_qiang1./2 ).^2 ./ 2 ./ (sigma_I_0L_qiang1));
       
% Pr_I_qiang2 = 1./(sqrt(2*pi.*sigma_I_0L_qiang2))./(I).*exp(-(log(I./I_0L) + ...
%            sigma_I_0L_qiang2./2+2*r.^2/W.^2).^2./(2.*sigma_I_0L_qiang2));
       
Pr_I_qiang2 = (2 .* pi .* sigma_I_0L_qiang2 .* I .* I).^(-0.5) .* ...
exp(-(log(I ./ I_0L) + sigma_I_0L_qiang2./2 ).^2 ./ 2 ./ (sigma_I_0L_qiang2));       
       
%弱湍流情况下
%概率密度函数，最后算BER做积分   
% 
% Pr_I_ruo = 1./(sqrt(2*pi.*sigma_I_0L_ruo))./(I).*exp(-(log(I./I_0L) + ...
%            sigma_I_0L_ruo./2+2*r.^2/W.^2).^2./(2.*sigma_I_0L_ruo));
       
 Pr_I_ruo = (2 .* pi .* sigma_I_0L_ruo .* I .* I).^(-0.5) .* ...
exp(-(log(I ./ I_0L) + sigma_I_0L_ruo./2 ).^2 ./ 2 ./ (sigma_I_0L_ruo));        
       
%夜视模式的大气湍流情况


% Pr_I_dark = 1./(sqrt(2*pi.*sigma_I_0L_dark))./(I).*exp(-(log(I./I_0L) + ...
%            sigma_I_0L_dark./2+2*r.^2/W.^2).^2./(2.*sigma_I_0L_dark));
       
 Pr_I_dark = (2 .* pi .* sigma_I_0L_dark .* I .* I).^(-0.5) .* ...
exp(-(log(I ./ I_0L) + sigma_I_0L_dark./2 ).^2 ./ 2 ./ (sigma_I_0L_dark));
%验证
 
prI = @(I) 1./(sqrt(2*pi.*sigma_I_0L_dark))./(I).*exp(-(log(I./I_0L) + ...
             sigma_I_0L_dark./2).^2./(2.*sigma_I_0L_dark));
% prI = @(I) (2 .* pi .* sigma_I_0L_dark .* I .* I).^(-0.5) .* ...
% exp(-(log(I ./ I_0L) + sigma_I_0L_dark./2 ).^2 ./ 2 ./ (sigma_I_0L_dark).^2);
% prI = @(I) 1./(sqrt(2*pi.*sigma_I_0L_ruo))./(I).*exp(-(log(I./I_0L) + ...
%         sigma_I_0L_ruo./2+2*r.^2/W.^2).^2./(2.*sigma_I_0L_ruo));
PIji = integral(prI,10^-15,3*10^-5)


%画图调试

% figure;
% plot(I,Pr_I_qiang1);
% 
% hold on ; 
% plot(I,Pr_I_qiang2);
% plot(I,Pr_I_ruo);
% plot(I,Pr_I_dark);
% %plot(I,log10(Pr_I_ruo));                          %用来确定I的积分范围
% title("强湍流和弱湍流,夜视情况下的PDF");
% legend("强湍流","弱湍流","夜视");
% hold off;

% figure;
% plot(I,log10(Pr_I_qiang1));
% hold on ;
% plot(I,log10(Pr_I_qiang2));
% plot(I,log10(Pr_I_ruo));       
% plot(I,log10(Pr_I_dark));  %用来确定I的积分范围
% title("强湍流和弱湍流,夜视情况下的PDF(log)");
% legend("强湍流","弱湍流","夜视");
% hold off;


%%
%------------------------------混沌失配计算与方差推导推导------------------------------%


yita1 = 8;%发射端PD探测效率
yita2 = 8;%发射端PD探测效率
RF1 = 100;%发射端RF放大倍数
RF2 = 100;%接收端RF放大倍数
V1 = 4.2;%发射端VpiRF
V2 = 4.2;%接收端VpiRF

%beta 论文里的环路增益差值  直接给定值
% beta=10;                            %chaos feedback gain
% beta_1=9.54;                        %失配误差直接给值
% deltabeta=beta_1-beta;              %论文里的环路增益差值
GEDFA2 = 1/(GEDFA1*alphabeta*(P1/P2));
A1 = 100;
A2 = A1;
G_shuaijian = 10^-1.62;  %光衰减器      %-14db大小的衰减器
betaemi = pi*1/4*RF1*yita1*P1*A1/(2*V1)*G_shuaijian;
betarec = pi*1/2*RF2*yita2*A2*GEDFA2./(1+maskeff) .* I ./ (2*V2)*G_shuaijian;
deltabetarate = abs(betaemi - betarec) ./ betaemi;

%论文种的K电流大小

K=I *R_d *GEDFA2*G_shuaijian/(1+maskeff)*G^2;           %*Ts               %论文里的大气端K
K_1=R_d*P2*G_shuaijian*G^2;                         %论文里的LD2的电流K
deltaK=K_1 - K;                       %论文里的电路(光功率)差值

deltaTdelay=0;                     %失配时延
tao=25*10^(12);                    %高截止响应时间
deltatao=0.25*10^(12);             %高不匹配中断响应时间。


% deltaK=0.002; 
deltafai=0.02;
%deltaK/K_1
%synchronization error 同步误差
yipsen_2=1/3*(deltaTdelay./tao).^2+(deltabetarate).^2+(1-pi/4)*(deltatao./tao).^2....
-2*(1-pi/4)*deltabetarate.*deltatao./tao-2*(1-pi/4)*deltaTdelay./tao*deltatao./tao;
%deltaK./K_1
n_rms=1/2.*K.^2.*(yipsen_2+deltafai^2+1/4*(deltaK/K_1).^2);  %失配噪声

% yipsen = 1/3*(deltaTdelay./tao)^2+(deltabetarate).^2+(1-pi/4)*(deltatao./tao)^2....
% -2*(1-pi/4)*deltabetarate.*deltatao./tao-2*(1-pi/4)*deltaTdelay./tao*deltatao./tao;%同步误差方差 与K无关
% deltamm = 0.5 * (K .^ 2) .* (yipsen + deltafai^2+0.25*(deltaK/K_1).^2);%失配方差


AAA=(yipsen_2 + deltafai^2+0.25*(deltaK/K_1).^2);
 %梯形积分，计算相关BER
 
pr_I  = Pr_I_dark';
pr_I1 = Pr_I_ruo';
pr_I2 = Pr_I_qiang1';
pr_I3 = Pr_I_qiang2';
% BPSK_BER_qiang2(1,i)=BER(i,:)*[pr_I3(1)./2;pr_I3(2:num-1);pr_I3(num)./2].* dI;  %梯形积分，计算相关BER

  hundun_dark = n_rms(1,:)*[pr_I(1)./2;pr_I(2:num-1);pr_I(num)./2].* dI;
  hundun_ruo = n_rms(1,:)*[pr_I1(1)./2;pr_I1(2:num-1);pr_I1(num)./2].* dI;
  hundun_qiang1 = n_rms(1,:)*[pr_I2(1)./2;pr_I2(2:num-1);pr_I2(num)./2].* dI;
  hundun_qiang2 = n_rms(1,:)*[pr_I3(1)./2;pr_I3(2:num-1);pr_I3(num)./2].* dI;
  
  loghundun_dark = log10(hundun_dark);
  loghundun_ruo= log10(hundun_ruo);
  loghundun_qiang1= log10(hundun_qiang1);
  loghundun_qiang2= log10(hundun_qiang2);
%plot(I,n_rms);                      %不是以P2为中心 会产生一定的偏移
%%
% plot(I,AAA);   
% figure;
% plot(I,Pr_I_dark); 
%%
%-------------------------------  混沌失配参数的计算   -------------------------------%
%I是接收光的分布
% P_chao_mismatch = abs(I - P2);                          %失配后混沌相消的光功率 1*100001的矩阵
%P_chao_mismatch = P2 - I;    
%plot(I,P_chao_mismatch)   debug用的 没用
P_chao_mismatch = I * GEDFA2 *maskeff./(1+maskeff)*G_shuaijian;       %jishou光强
P_jieshou = P_chao_mismatch*[pr_I(1)./2;pr_I(2:num-1);pr_I(num)./2].* dI;
%%

%%
%------------------------------BPSK调制后的信号功率计算--------------------------------%
%探测到的信号大小
signal_P1 = (1/2*(P2+P_chao_mismatch)+sqrt(P2*P_chao_mismatch)) ; %I1的功率
signal_P2 = (1/2*(P2+P_chao_mismatch)-sqrt(P2*P_chao_mismatch)) ; %I2的功率

signal_power = 2*G^2*R_d^2*P_chao_mismatch*P2; % debug用的数，没有用
%%
signal_power_jieshou = signal_power*[pr_I(1)./2;pr_I(2:num-1);pr_I(num)./2].* dI;
%%
%------------------噪声功率计算 探测器噪声计算 热噪声和散粒噪声--------------------%
% 散粒噪声和热噪声
% 散粒噪声
sita_s1=2*G^2*e*F_A*delt_F*(R_d.*signal_P1+I_d);  %上面的散粒噪声
sita_s2=2*G^2*e*F_A*delt_F*(R_d.*signal_P2+I_d);  %下面的散粒噪声
%(1/2*(P_LO+PT_h)+sqrt(P_LO*PT_h)    是信号I

% 热噪声
sita_T=4*k_B*T*F_n*delt_F/R_L;

%相位噪声
phi = -0.5:0.001:0.5;                            %激光器相位噪声



%%
% 
% %--------------------------------静态BER的计算-----------------------------------%

delta_v=linspace(1,1200,400)*10^3;  
BER_p = zeros(length(signal_power),length(delta_v));
beiji = zeros(1,length(phi));

bc1 = (max(phi)-min(phi))/(length(phi)-1);
phase_noise= zeros(1,length(delta_v));
noise =zeros(length(signal_power),length(delta_v));
SNR = zeros(length(signal_power),length(delta_v));   
for j = 1:length(signal_power)
phase_noise = signal_power(j)./2.*(1-exp(-8*pi*delta_v*Ts));
noise(j,:) =(sita_s1(j)+sita_s2(j)+2*sita_T+phase_noise+n_rms(j));%
SNR(j,:) = signal_power(j)./(sita_s1(j)+sita_s2(j)+2*sita_T+phase_noise+n_rms(j));%
BER_p(j,:) =1/2*erfc((SNR(j,:)/2).^0.5);    
end

%%
phase_noise1= zeros(length(delta_v),length(signal_power));
for j = 1:length(signal_power)
phase_noise1(:,j) = signal_power(j)./2.*(1-exp(-8*pi*delta_v*Ts));
end

%%
phase_dark =zeros(1,length(delta_v));
phase_ruo =zeros(1,length(delta_v));
phase_qiang1 =zeros(1,length(delta_v));
phase_qiang2 =zeros(1,length(delta_v));
for i =1:length(delta_v)
phase_dark(1,i) = phase_noise1(i,:)*[pr_I(1)./2;pr_I(2:num-1);pr_I(num)./2].* dI;
phase_ruo(1,i) = phase_noise1(i,:)*[pr_I1(1)./2;pr_I1(2:num-1);pr_I1(num)./2].* dI;
phase_qiang1(1,i) = phase_noise1(i,:)*[pr_I2(1)./2;pr_I2(2:num-1);pr_I2(num)./2].* dI;
phase_qiang2 (1,i)= phase_noise1(i,:)*[pr_I3(1)./2;pr_I3(2:num-1);pr_I3(num)./2].* dI;
end

%%
%计算比值的部分

rate_dark =  hundun_dark ./phase_dark;
rate_ruo =  hundun_ruo ./ phase_ruo;
rate_qiang1 =  hundun_qiang1 ./ phase_qiang1;
rate_qiang2 =  hundun_qiang2 ./phase_qiang2;

rate_fianl = [rate_dark;rate_ruo;rate_qiang1;rate_qiang2];


figure ;
for i = 1:4
subplot(2,2,i)
plot(delta_v/1000,rate_fianl(i,:),'c','LineWidth',2.5);
set(gca,'linewidth',2.5)
xlabel('Linewidth \Deltav/kHz');
ylabel('Noise ratio');
set(gca,'FontSize',18);
ylim([0 7]);
set(gca,'YTick',[0:1:7]);
if(i == 1)
    legend('\sigma_I^2=0.05','FontSize',18');
elseif i == 2
     legend('\sigma_I^2=0.10','FontSize',18');    
  
elseif i == 3
     legend('\sigma_I^2=0.15','FontSize',18');    
 
elseif i == 4
     legend('\sigma_I^2=0.20','FontSize',18');    
end
set(gca,'Fontname','times new Roman');
set(gca,'xtick',[0 200 400 600 800 1000]);

end
%%
% figure;
% hold on;
% plot(delta_v,rate_dark);
% plot(delta_v,rate_ruo);
% plot(delta_v,rate_qiang1);
% plot(delta_v,rate_qiang2);
% hold off;
% ylim([0 8]);
% 
% %%
% figure ;
% for i = 1:4
% subplot(2,2,i)
% plot(delta_v,log10(rate_dark));
% plot(delta_v,log10(rate_ruo));
% plot(delta_v,log10(rate_qiang1));
% plot(delta_v,log10(rate_qiang2));
% end
% 
% %%
% %%画图
% figure;
% b=bar3(rate_fianl,0.6);
% %  set(gca,'ZDir','reverse') ;
% colorbar;
% %%
% % caxis([-9,-4]);
% 
% % zlim([-9,-4]);      %z坐标轴范围
% for k = 1:length(b)
%     zdata = b(k).ZData;
%     b(k).CData = zdata;
%     b(k).FaceColor = 'interp';
% end
% 
% xlabel('Linewidth \Delta\upsilon / kHz','FontWeight','bold','Rotation',28,'Position',[10 3 0],'FontName','Times New Roman');
% ylabel('\sigma_I^2 ','Rotation',-32,'Position',[2 -0.5 0],'FontName','Times New Roman');
% zlabel('phase noise power ','Rotation',-90,'Position',[-1 0 -8]);
% set(gca,'FontSize',22,'FontWeight','bold');
% set(gca,'linewidth',3);
% set(gca,'xticklabel',{'1','60','120','180','240','300','360','420','480','540','600','660','720','780','900'});
% set(gca,'yticklabel',{'0.05','0.10','0.15','0.20'});
% set(gca,'XGrid', 'off', 'YGrid', 'off','ZGrid', 'off');
% set(gca,'Fontname','times new Roman');

% ylim([0,4.5]);  


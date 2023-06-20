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

PIji = integral(prI,10^-15,3*10^-5)




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
deltaK=0;                       %论文里的电路(光功率)差值

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

%%

figure;
box on;
plot(I(1,1:450),deltabetarate(1,1:450),'-','LineWidth',4,'Color',[0 0 0.6]);
xlabel('optical power /W','FontWeight','bold','FontName','Times New Roman','FontSize',35);
ylabel('\Delta\beta /\beta ','FontWeight','bold','FontName','Times New Roman','FontSize',35);
set(gca,'FontSize',35,'FontWeight','bold');
set(gca,'linewidth',4);
set(gca,'Fontname','times new Roman');
%%

%画图
figure;
box on;
yyaxis left;
hold on ;
plot(I(1,1:450),Pr_I_dark(1,1:450),'-','LineWidth',4,'Color',[0.70 0.55 0.4]);
plot(I(1,1:450),Pr_I_ruo(1,1:450),'-','LineWidth',4,'Color',[0.2 0.6 0.1]);
plot(I(1,1:450),Pr_I_qiang1(1,1:450),'-','LineWidth',4,'Color',[0.02 0.7 0.7]);
plot(I(1,1:450),Pr_I_qiang2(1,1:450),'-','LineWidth',4,'Color',[0.8 0.02 0.8]);
xlabel('optical power /W','FontWeight','bold','FontName','Times New Roman','linewidth',4,'FontSize',35);
ylabel('PDF','FontWeight','bold','FontName','Times New Roman','linewidth',4,'FontSize',35);
set(gca,'FontSize',35,'FontWeight','bold');
set(gca,'linewidth',4);
set(gca,'Fontname','times new Roman');




 yyaxis right;
 
plot(I,n_rms,'r--','LineWidth',4);

 ylim([10^-12 10^-9]);
xlabel('optical power I /W','FontWeight','bold','FontName','Times New Roman','FontSize',35);
ylabel('chaotic noise power /W','FontWeight','bold','FontName','Times New Roman','FontSize',35);
set(gca,'FontSize',35,'FontWeight','bold');
set(gca,'linewidth',4);
set(gca,'Fontname','times new Roman');


legend('\sigma_I^2=0.05','\sigma_I^2=0.10','\sigma_I^2=0.15',...
'\sigma_I^2=0.20','chaotic noise power','NumColumns',1,'FontSize',35');



hold off;

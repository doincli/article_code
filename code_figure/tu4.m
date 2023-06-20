%---------------------------------������Ϣ---------------------------------%
% �������ƣ�bpsk_ber_v6
% ������Ϣ��honney_li
% ʵ�ֹ��ܣ����������������硢��λ��������µ�BPSK������
% ���뻷����matlab2018b
% ��д���ڣ�2022/10/12
%-------------------------------------------------------------------------%

%-------------------------------�汾˵��-------------------------------%

%{

1.��V2�汾�Ļ����ϼ������������BER   
2.������4�ּ��������߿��ֱ�Ϊ500 800 1200 2000KHz
3.��V3�汾�Ļ����ϼ��� ǿ���� ������ ����ҹ������µĴ������� PDF ���ػ��� ����һ��pdf������0
4.��V4�汾�Ļ����ϼ�����������Ĺ�ʽ����
5.����汾˵����word ������ַ��
https://artic-1302263000.cos.ap-nanjing.myqcloud.com/BER%E6%8E%A8%E5%AF%BC%E6%94%B93.docx
6.����������ع�ʽ�����ݹ�ʽ�޸ĳ���

%}

%-------------------------------�汾����-------------------------------%
%{

1.beta��ȡֵ���⣬��ѧ����ʽ�Ƕ�ֵ���������Ĵ����ǻ�����ʽ ��ѧ��С���Ĺ�ʽ11   
2.SNR��Ӱ�����⣬�����ĺ����͹�����ܶȺ�����ƥ��
3.��ʧ��������ʱ��û���õ���Ĥ�ȣ������������������Ҫ����
4.Ҫ�������༤�������߿����߿�Ҳ�ֲ����֣�����ֻ����4��

%}

%-------------------------------��������-------------------------------%
clear all;
close all;
clc;
format long;%��Ч����16λ

%-------------------------------ϵͳ����-------------------------------%
%һЩ������
hp = 6.6260693 * 10^(-34);          %���ʿ˳��� 
e = 1.60217733 * 10^(-19);          %���ӵ���
kc = 1.3806505 * 10^(-23);          %������������
c = 3 * 10^8;                       %����

%-------------------------------��·��ز���-------------------------------%
zenith_angle = 0;                   %�춥�ǣ���λ����
zeta = deg2rad(zenith_angle);       %���춥�ǣ���λ��rad
H = 38.*10.^6;                      %�����Ǿ������ĸ߶ȣ���λ:m
h0 = 100;                           %���������ǽ���վ�߶ȣ����������Ƿ���˸߶ȣ���λ��m
L=(H-h0).*sec(zeta);                %��·���룬��λ��m
alpha = 1;                          %����·���������ϵ����Ϊ1����û�����

%-------------------------------����˲���-------------------------------%
P1 = 8*1e-3;                       %���Ͷ˼������㶨����
lambda = 1550.*10^(-9);             %�˲�������λ��m
k = 2 .* pi ./ lambda;              %��������λ��rad/m
v = c ./ lambda;                    %��= c/�ˣ���Ƶ�ʣ���λ��Hz
%divergence_angle = 1;%��ɢ�ǣ�ȫ�ǣ�����λ����
%theta = deg2rad(divergence_angle);%����ɢ�ǣ�ȫ�ǣ�����λ��rad
theta = 30  .* 10^(-6);             %����ɢ�ǣ�ȫ�ǣ�����λ��rad
W0 = 0.1;                           %����ھ�,�����ߴ�С����λ��m����֤��߰뾶���������뾶
GEDFA1 = 2000;                      %����˵�EDFA����
maskeff = 1.2;                      %��Ĥ�� �źŹ��ʱ����ز�����

%-------------------------------���ƶ˲���-------------------------------%
Ps = maskeff * P1 ;                       %���ƶ� �źż������㶨���� 

%-------------------------------���ն˲���-------------------------------%
P2 = P1/2;                        %���ն˼������㶨����,���书�ʵ�һ�룬Ϊ�˻�������
W = W0 + theta .* L ./ 2;          %�������ϵĹ�߰뾶 W = ��L/2��Ϊ���ƹ�ʽ
Dr = 0.40;                         %���տھ�,��λ��m
Waist=2 .* lambda ./ (pi.*theta);  %�����뾶 
f = pi * (Waist.^2) / lambda;      %��������
Rr = L + f .^ 2 ./ L;              %�������ϵĲ�ǰ���ʰ뾶
%GEDFA2 = 400;                     %���ն˵�EDFA����

                                     
%-------------------------------BPSK���ģ����ز���-------------------------------%

%-------------------------------ʧ����������-------------------------------%

%beta=10;                           %chaos feedback gain ��·����

%-------------------------------̽��������-------------------------------%
%%ƽ��̽��������
yita=0.8;                        %������Ч��
R_d=yita*e/(hp*v);                %����PD��Ӧ��R_d

% ̽����������������
G=1;                             %̽��������
F_A=G^0.5;
I_d=10^-9;                       %������
delt_F=2*10^9;                 %������Ч����
Ts=1/delt_F;                     %������Чʱ��
T=300;                           %̽�����¶�
R_L=50;                          %̽�������ص���
k_B=1.380649*10^-23;             %������������
F_n=1;                           %̽������������


%-------------------------------  ������·��ǿ����  -------------------------------%
%r=0���չ�����ƽ�����ʴ�С ����������·��ģ����ն����Ĺ⹦��
Pt = (P1 + Ps);                                    %�����ڸǱȵĹ⹦��
r = 0;                                               %�����Ƕ�׼���
I_0L = alpha .* Pt .*GEDFA1.* Dr.^2 ./ (2 .* W.^2);  %������·�����,����EDFA2֮ǰ
 
alphabeta = Dr.^2 ./ (2 .* W.^2);

%-------------------------------  ������ǿ��˸��������   -------------------------------%


sigma_I_0L_dark =0.05;
sigma_I_0L_ruo=0.10; %sigma_I_0L������ǿ��˸���� ������

sigma_I_0L_qiang1=0.15; %sigma_I_0L������ǿ��˸���� ǿ����
sigma_I_0L_qiang2=0.20;
%����Ӱ���µĹ�ǿ�����ܶȺ���
num = 1000+1;

I = linspace(10^-15,3*10^-5,num);               %���չ����幦�ʴ�С,Ϊ����  Ϊ������Ȩ����BER
dI =(max(I)-min(I))./num;                      %����

%-------------------------------  ������·�����ܶȺ�������   -------------------------------%
%ǿ���������

Pr_I_qiang_test = 1./(sqrt(2*pi.*sigma_I_0L_qiang1))./(I).*exp(-(log(I./I_0L) + ...
            sigma_I_0L_qiang1./2).^2./(2.*sigma_I_0L_qiang1));
       
Pr_I_qiang1 = (2 .* pi .* sigma_I_0L_qiang1 .* I .* I).^(-0.5) .* ...
                 exp(-(log(I ./ I_0L) + sigma_I_0L_qiang1./2 ).^2 ./ 2 ./ (sigma_I_0L_qiang1));
       
% Pr_I_qiang2 = 1./(sqrt(2*pi.*sigma_I_0L_qiang2))./(I).*exp(-(log(I./I_0L) + ...
%            sigma_I_0L_qiang2./2+2*r.^2/W.^2).^2./(2.*sigma_I_0L_qiang2));
       
Pr_I_qiang2 = (2 .* pi .* sigma_I_0L_qiang2 .* I .* I).^(-0.5) .* ...
exp(-(log(I ./ I_0L) + sigma_I_0L_qiang2./2 ).^2 ./ 2 ./ (sigma_I_0L_qiang2));       
       
%�����������
%�����ܶȺ����������BER������   
% 
% Pr_I_ruo = 1./(sqrt(2*pi.*sigma_I_0L_ruo))./(I).*exp(-(log(I./I_0L) + ...
%            sigma_I_0L_ruo./2+2*r.^2/W.^2).^2./(2.*sigma_I_0L_ruo));
       
 Pr_I_ruo = (2 .* pi .* sigma_I_0L_ruo .* I .* I).^(-0.5) .* ...
exp(-(log(I ./ I_0L) + sigma_I_0L_ruo./2 ).^2 ./ 2 ./ (sigma_I_0L_ruo));        
       
%ҹ��ģʽ�Ĵ����������


% Pr_I_dark = 1./(sqrt(2*pi.*sigma_I_0L_dark))./(I).*exp(-(log(I./I_0L) + ...
%            sigma_I_0L_dark./2+2*r.^2/W.^2).^2./(2.*sigma_I_0L_dark));
       
 Pr_I_dark = (2 .* pi .* sigma_I_0L_dark .* I .* I).^(-0.5) .* ...
exp(-(log(I ./ I_0L) + sigma_I_0L_dark./2 ).^2 ./ 2 ./ (sigma_I_0L_dark));
%��֤
 
prI = @(I) 1./(sqrt(2*pi.*sigma_I_0L_dark))./(I).*exp(-(log(I./I_0L) + ...
             sigma_I_0L_dark./2).^2./(2.*sigma_I_0L_dark));

PIji = integral(prI,10^-15,3*10^-5)




%%
%------------------------------����ʧ������뷽���Ƶ��Ƶ�------------------------------%


yita1 = 8;%�����PD̽��Ч��
yita2 = 8;%�����PD̽��Ч��
RF1 = 100;%�����RF�Ŵ���
RF2 = 100;%���ն�RF�Ŵ���
V1 = 4.2;%�����VpiRF
V2 = 4.2;%���ն�VpiRF

%beta ������Ļ�·�����ֵ  ֱ�Ӹ���ֵ
% beta=10;                            %chaos feedback gain
% beta_1=9.54;                        %ʧ�����ֱ�Ӹ�ֵ
% deltabeta=beta_1-beta;              %������Ļ�·�����ֵ
GEDFA2 = 1/(GEDFA1*alphabeta*(P1/P2));
A1 = 100;
A2 = A1;
G_shuaijian = 10^-1.62;  %��˥����      %-14db��С��˥����
betaemi = pi*1/4*RF1*yita1*P1*A1/(2*V1)*G_shuaijian;
betarec = pi*1/2*RF2*yita2*A2*GEDFA2./(1+maskeff) .* I ./ (2*V2)*G_shuaijian;
 deltabetarate = abs(betaemi - betarec) ./ betaemi;

%�����ֵ�K������С

K=I *R_d *GEDFA2*G_shuaijian/(1+maskeff)*G^2;           %*Ts               %������Ĵ�����K
K_1=R_d*P2*G_shuaijian*G^2;                         %�������LD2�ĵ���K
deltaK=0;                       %������ĵ�·(�⹦��)��ֵ

deltaTdelay=0;                     %ʧ��ʱ��
tao=25*10^(12);                    %�߽�ֹ��Ӧʱ��
deltatao=0.25*10^(12);             %�߲�ƥ���ж���Ӧʱ�䡣


% deltaK=0.002; 
deltafai=0.02;
%deltaK/K_1
%synchronization error ͬ�����
yipsen_2=1/3*(deltaTdelay./tao).^2+(deltabetarate).^2+(1-pi/4)*(deltatao./tao).^2....
-2*(1-pi/4)*deltabetarate.*deltatao./tao-2*(1-pi/4)*deltaTdelay./tao*deltatao./tao;
%deltaK./K_1
n_rms=1/2.*K.^2.*(yipsen_2+deltafai^2+1/4*(deltaK/K_1).^2);  %ʧ������

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

%��ͼ
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

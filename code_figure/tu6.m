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
% prI = @(I) (2 .* pi .* sigma_I_0L_dark .* I .* I).^(-0.5) .* ...
% exp(-(log(I ./ I_0L) + sigma_I_0L_dark./2 ).^2 ./ 2 ./ (sigma_I_0L_dark).^2);
% prI = @(I) 1./(sqrt(2*pi.*sigma_I_0L_ruo))./(I).*exp(-(log(I./I_0L) + ...
%         sigma_I_0L_ruo./2+2*r.^2/W.^2).^2./(2.*sigma_I_0L_ruo));
PIji = integral(prI,10^-15,3*10^-5)


%��ͼ����

% figure;
% plot(I,Pr_I_qiang1);
% 
% hold on ; 
% plot(I,Pr_I_qiang2);
% plot(I,Pr_I_ruo);
% plot(I,Pr_I_dark);
% %plot(I,log10(Pr_I_ruo));                          %����ȷ��I�Ļ��ַ�Χ
% title("ǿ������������,ҹ������µ�PDF");
% legend("ǿ����","������","ҹ��");
% hold off;

% figure;
% plot(I,log10(Pr_I_qiang1));
% hold on ;
% plot(I,log10(Pr_I_qiang2));
% plot(I,log10(Pr_I_ruo));       
% plot(I,log10(Pr_I_dark));  %����ȷ��I�Ļ��ַ�Χ
% title("ǿ������������,ҹ������µ�PDF(log)");
% legend("ǿ����","������","ҹ��");
% hold off;


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
deltaK=K_1 - K;                       %������ĵ�·(�⹦��)��ֵ

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

% yipsen = 1/3*(deltaTdelay./tao)^2+(deltabetarate).^2+(1-pi/4)*(deltatao./tao)^2....
% -2*(1-pi/4)*deltabetarate.*deltatao./tao-2*(1-pi/4)*deltaTdelay./tao*deltatao./tao;%ͬ������ ��K�޹�
% deltamm = 0.5 * (K .^ 2) .* (yipsen + deltafai^2+0.25*(deltaK/K_1).^2);%ʧ�䷽��


AAA=(yipsen_2 + deltafai^2+0.25*(deltaK/K_1).^2);
 %���λ��֣��������BER
 
pr_I  = Pr_I_dark';
pr_I1 = Pr_I_ruo';
pr_I2 = Pr_I_qiang1';
pr_I3 = Pr_I_qiang2';
% BPSK_BER_qiang2(1,i)=BER(i,:)*[pr_I3(1)./2;pr_I3(2:num-1);pr_I3(num)./2].* dI;  %���λ��֣��������BER

  hundun_dark = n_rms(1,:)*[pr_I(1)./2;pr_I(2:num-1);pr_I(num)./2].* dI;
  hundun_ruo = n_rms(1,:)*[pr_I1(1)./2;pr_I1(2:num-1);pr_I1(num)./2].* dI;
  hundun_qiang1 = n_rms(1,:)*[pr_I2(1)./2;pr_I2(2:num-1);pr_I2(num)./2].* dI;
  hundun_qiang2 = n_rms(1,:)*[pr_I3(1)./2;pr_I3(2:num-1);pr_I3(num)./2].* dI;
  
  loghundun_dark = log10(hundun_dark);
  loghundun_ruo= log10(hundun_ruo);
  loghundun_qiang1= log10(hundun_qiang1);
  loghundun_qiang2= log10(hundun_qiang2);
%plot(I,n_rms);                      %������P2Ϊ���� �����һ����ƫ��
%%
% plot(I,AAA);   
% figure;
% plot(I,Pr_I_dark); 
%%
%-------------------------------  ����ʧ������ļ���   -------------------------------%
%I�ǽ��չ�ķֲ�
% P_chao_mismatch = abs(I - P2);                          %ʧ�����������Ĺ⹦�� 1*100001�ľ���
%P_chao_mismatch = P2 - I;    
%plot(I,P_chao_mismatch)   debug�õ� û��
P_chao_mismatch = I * GEDFA2 *maskeff./(1+maskeff)*G_shuaijian;       %jishou��ǿ
P_jieshou = P_chao_mismatch*[pr_I(1)./2;pr_I(2:num-1);pr_I(num)./2].* dI;
%%

%%
%------------------------------BPSK���ƺ���źŹ��ʼ���--------------------------------%
%̽�⵽���źŴ�С
signal_P1 = (1/2*(P2+P_chao_mismatch)+sqrt(P2*P_chao_mismatch)) ; %I1�Ĺ���
signal_P2 = (1/2*(P2+P_chao_mismatch)-sqrt(P2*P_chao_mismatch)) ; %I2�Ĺ���

signal_power = 2*G^2*R_d^2*P_chao_mismatch*P2; % debug�õ�����û����
%%
signal_power_jieshou = signal_power*[pr_I(1)./2;pr_I(2:num-1);pr_I(num)./2].* dI;
%%
%------------------�������ʼ��� ̽������������ ��������ɢ������--------------------%
% ɢ��������������
% ɢ������
sita_s1=2*G^2*e*F_A*delt_F*(R_d.*signal_P1+I_d);  %�����ɢ������
sita_s2=2*G^2*e*F_A*delt_F*(R_d.*signal_P2+I_d);  %�����ɢ������
%(1/2*(P_LO+PT_h)+sqrt(P_LO*PT_h)    ���ź�I

% ������
sita_T=4*k_B*T*F_n*delt_F/R_L;

%��λ����
phi = -0.5:0.001:0.5;                            %��������λ����



%%
% 
% %--------------------------------��̬BER�ļ���-----------------------------------%

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
%�����ֵ�Ĳ���

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
% %%��ͼ
% figure;
% b=bar3(rate_fianl,0.6);
% %  set(gca,'ZDir','reverse') ;
% colorbar;
% %%
% % caxis([-9,-4]);
% 
% % zlim([-9,-4]);      %z�����᷶Χ
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


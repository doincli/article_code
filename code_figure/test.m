% ��ָ������µ��߿�����
% С����ͼ2��ͼ3

tic
% close all
% clear all
%% ��������
lambda = 1.55*10^-6;
k = 2*pi/lambda;
zeta = 0;           %�춥��=0
H = 3.8*10^7;       %���Ǹ߶�Ϊ38000km
h0 = 100;           %������ն˸߶�
syms h;             %����߶Ȳ���h������Ӱ������ha����������ƫ�����r
ha = 0.1:0.1:20;
r = 1:1:1000;
% % phi=-3:1:3;

w = 21;             %rms����Ϊ21m/s
A0 = 1.7*10^-14;

% ̽������Ӧ��
yita=0.75;
e=1.6*10^-19;
hp=6.626*10^-34;     %���ʿ˳���
c=3*10^8;
v=c/lambda;
R_d=yita*e/(hp*v);   %����PD��Ӧ��R_d

% ��·����
a_loss=1;           %�������ϵ��
L=3.8*10^7;         %��·����38000km
D_r=0.6;            %���տ׾�0.6m
theta=30*10^-6;     %��ɢ��30΢rad
W0=0.1;             %���Ͷ˹����뾶0.1m
W=W0+theta*L/2;     %���ն˹����뾶
P1=1*10^-3;         %���Ͷ˼�����LD1�������Ϊ1mW
G_EDFA=150;         %EDFA�Ŵ���Ϊ1000        
a_link=a_loss*(D_r)^2/(2*W^2);   %��·�������ϵ��
% PT = 2*10^-8;       %δ����ָ����·��ĵĴ��½��չ���

% ̽����������������
G=1;
F_A=G^0.5;
I_d=10^-9;
delt_F=2.5*10^9;
Ts=1/delt_F;
T=300;
R_L=50;
k_B=1.380649*10^-23;
F_n=1;
P_LO=0.001;         %����⹦��

%% gamma-gamma����
% ��������ǿ��sigma2_R
Cn = 0.00594*((w/27)^2)*((10^-5*h)^10)*exp(-h/1000)+(2.7*10^-16)*exp(-h/1500)+A0*exp(-h/100);
x = Cn*(h-h0)^(5/6);
fun1 = matlabFunction(x);
y = integral(@(h)fun1(h),h0,H);
sigma2_R = 2.25*k^(7/6)*(sec(zeta))^(11/6)*2*y;         %Rytov�����ʾ����ǿ��

% ����ǿ�������������߶�
a = (exp(0.49*sigma2_R/(1+1.1*sigma2_R^(12/5))^(7/6))-1)^-1;     %aΪ�������߶�
b = (exp(0.51*sigma2_R/(1+0.69*sigma2_R^(12/5))^(5/6))-1)^-1;    %bΪǿ�����߶�

% ha��PDF
be1 = 2*(a*b*ha).^(1/2);      %�ڶ���������������������Ĳ���
fh = (2*(a*b)^((a+b)/2)/(gamma(a)*gamma(b)))*(ha).^((a+b)/2-1).*besselk(a-b,be1); %ha��PDF
% figures
% plot(ha,fh);

%% ������·�������
PT=G_EDFA*a_link*P1;        %����������·��ģ����ն����Ĺ⹦��
% PT_r=PT*exp(-2*r.^2/W^2);   %���ն˹⹦����r�ķֲ�


%% ���������Խ��չ⹦��P��Ӱ��
bc = (max(ha)-min(ha))/(length(ha)-1);
fh = fh';
PT_h = PT*ha*[fh(1)./2;fh(2:length(ha)-1);fh(length(ha))./2].*bc;

% beiji = PT*ha.*fh;
% bc = (max(ha)-min(ha))/(length(ha)-1);
% pp = beiji.*bc;
% PT_h = sum(pp(:));
I_ac=2*G^2*R_d^2*PT_h*P_LO;

%% ɢ��������������
% ɢ������
sita_s1=2*G^2*e*F_A*delt_F*(R_d*(1/2*(P_LO+PT_h)+sqrt(P_LO*PT_h))+I_d);
sita_s2=2*G^2*e*F_A*delt_F*(R_d*(1/2*(P_LO+PT_h)-sqrt(P_LO*PT_h))+I_d);

% ������
sita_T=4*k_B*T*F_n*delt_F/R_L;

%% ��λ����
phi = -0.5:0.001:0.5;
% phi = single(phi);
% phi = gpuArray(phi);
delta_v=linspace(1,1200,200)*10^3;                 %�������߿�
BER = zeros(1,length(delta_v));

% ��ͬ�߿��µ���λ�����ֲ�
linewidth = [100 300 500 1000]*10^3;
fp1 = zeros(length(linewidth),length(phi));
for i = 1:length(linewidth)
sigma1 = 4*pi*linewidth(i)*Ts;
fp1(i,:) = (1/(sqrt(2*pi)*sqrt(sigma1)))*gaussmf(phi,[sqrt(sigma1) 0]);
end

% ����BER��SNR
bc1 = (max(phi)-min(phi))/(length(phi)-1);
for i = 1:length(delta_v)
    sigma_2=4*pi*delta_v(i)*Ts;                             %��λ��������
    fp=(1/(sqrt(2*pi)*sqrt(sigma_2)))*gaussmf(phi,[sqrt(sigma_2) 0]); %��λ��������PDF
    sita_PN = I_ac*(sin(phi)).^2;
    SNR=(I_ac./(sita_s1+sita_s2+2*sita_T+sita_PN));
    BER_p =1/2*erfc((SNR/2).^0.5);
    beiji=BER_p.*fp*bc1;
    BER(i) = sum(beiji(:));
end

for i = 1:length(linewidth)
    beiji2(i,:)= (BER_p.*fp1(i,:));
end

%% ��ͬ�߿���ber(theta)*f(theta,delta v )��ͼ
figure
%������ɫ����
fig1_color = [[0.8500 0.3250 0.0980]
              [0.4660 0.6740 0.1880]
              [0.9290 0.6940 0.1250]
              [0.4940 0.1840 0.5560]];
%����ͼ��
str = ['\Deltav = 100 kHz',...
        '\Deltav = 300 kHz',...
        '\Deltav = 500 kHz',...
        '\Deltav = 1 MHz'];
for i = 1:4
    subplot(2,2,i)
    hold on 
    plot(phi,beiji2(i,:),'Color',fig1_color(i,:),'LineStyle','-','linewidth',2);
    area(phi,beiji2(i,:),'FaceColor',fig1_color(i,:),'FaceAlpha',0.5,'edgecolor','none');
    box on
    xlabel('Phase noise \theta/rad');
    ylabel('\itBER(\theta)��f(\theta,\Deltav)');
    hold off
    set(gca,'LineWidth',1);
    set(gca,'Fontname','times new Roman');
    set(gca,'Fontweight','bold');
    set(gca,'Fontsize',10);
end

% fill([phi phi(1001) phi(1)],[beiji2(4,:) -120 -120],[0.4940 0.1840 0.5560],'facealpha',0.5,'edgecolor','none');
% fill([phi phi(1001) phi(1)],[beiji2(3,:) -120 -120],[0.9290 0.6940 0.1250],'facealpha',0.5,'edgecolor','none');
% fill([phi phi(1001) phi(1)],[beiji2(2,:) -120 -120],[0.4660 0.6740 0.1880],'facealpha',0.5,'edgecolor','none');
% fill([phi phi(1001) phi(1)],[beiji2(1,:) -120 -120],[0.8500 0.3250 0.0980],'facealpha',0.5,'edgecolor','none');


% legend({str(1),...
%         '\Deltav = 300 kHz',...
%         '\Deltav = 500 kHz',...
%         '\Deltav = 1 MHz',...
%          },'Fontsize',12,'Location','southeast');
% legend('boxoff');

%% ��ͬ�߿�����λ��������ͼ
P_LPN = sita_PN*50*1000;        %˲ʱ��λ��������
P_p = zeros(1,length(linewidth)); %������λ�����⹦��
for i = 1:length(linewidth)
    temp = P_LPN.*fp1(i,:)*bc1;
    P_p(i)= sum(temp(:));
end

figure
hold on
yyaxis right
plot(phi,10*log10(P_LPN),'Color',[0 0.4470 0.7410],'linewidth',2)
set(gca,'Ycolor','k');
ylabel('P_{LPN}(\theta)/dBm');

yyaxis left
for i = 1:4
    plot(phi,fp1(i,:),'Color',fig1_color(i,:),'LineStyle','-','linewidth',2)
end
set(gca,'Ycolor','k');
ylabel('PDF of \theta');
xlabel('Phase noise \theta/rad');
hold off
box on 
set(gca,'LineWidth',1);
set(gca,'Fontname','times new Roman');
set(gca,'Fontweight','bold');
set(gca,'Fontsize',12);

legend({'\Deltav = 100 kHz',...
        '\Deltav = 300 kHz',...
        '\Deltav = 500 kHz',...
        '\Deltav = 1 MHz',...
         'P_{LPN}(\theta)',},'Fontsize',12,'Location','southeast');
legend('boxoff');


%% ������λ��������
for i = 1:length(linewidth)
    beiji3(i,:)= (sita_PN*50*1000.*fp1(i,:));
end
figure
%������ɫ����
fig1_color = [[0.8500 0.3250 0.0980]
              [0.4660 0.6740 0.1880]
              [0.9290 0.6940 0.1250]
              [0.4940 0.1840 0.5560]];
%����ͼ��
for i = 1:4
    subplot(2,2,i)
    hold on 
    plot(phi,beiji3(i,:),'Color',fig1_color(i,:),'LineStyle','-','linewidth',2);
    area(phi,beiji3(i,:),'FaceColor',fig1_color(i,:),'FaceAlpha',0.5,'edgecolor','none');
    box on
    xlabel('Phase noise \theta/rad');
    ylabel('\itP_{pn}(\theta)��f(\theta,\Deltav)/mW');
    hold off
    set(gca,'LineWidth',1);
    set(gca,'Fontname','times new Roman');
    set(gca,'Fontweight','bold');
    set(gca,'Fontsize',12);
end
% figure
% plot(phi,SNR);
%% ��ͼ
% ��ͬ�߿��µ���λ����PDF

% BER vs ��λ�����Ͳ�ͬ�߿��µ���λ����PDF
figure 
yyaxis right
f0 = plot(phi,log10(BER_p),'r','linewidth',2);
f0.Color = '#0072BD';
ylabel('lg[BER(\theta)]');
set(gca,'Ycolor','k');
yyaxis left
box on
hold on 
f1 = plot(phi,fp1(1,:),'-','linewidth',2);
f2 = plot(phi,fp1(2,:),'-','linewidth',2);
f3 = plot(phi,fp1(3,:),'-','linewidth',2);
f4 = plot(phi,fp1(4,:),'-','linewidth',2);
f1.Color = '#D95319';
f2.Color = '#77AC30';
f3.Color = '#EDB120';
f4.Color = '#7E2F8E';
hold off 
ylabel('PDF of \theta');
xlabel('Phase noise \theta/rad');
legend({'\Deltav = 100 kHz',...
        '\Deltav = 300 kHz',...
        '\Deltav = 500 kHz',...
        '\Deltav = 1 MHz',...
        'lg(BER(\theta))',...
         },'Fontsize',12,'Location','southeast');
legend('boxoff');
set(gca,'Ycolor','k');


set(gca,'LineWidth',1);
set(gca,'Fontname','times new Roman');
set(gca,'Fontweight','bold');
set(gca,'Fontsize',12);

% BER vs�߿�ͼ
figure
% subplot(1,2,2)
plot(delta_v/1000,log10(BER),'linewidth',2);
xlabel('Linewidth \Delta{\itv}/kHz');
ylabel('lg(BER)');

set(gca,'LineWidth',1);
% set(gca,'YTick',-15:5:-5);
set(gca,'Fontname','times new Roman');
set(gca,'Fontweight','bold');
set(gca,'Fontsize',12);



toc

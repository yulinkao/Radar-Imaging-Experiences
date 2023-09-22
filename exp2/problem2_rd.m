%名?：基于RD算法的雷?成像仿真

close all;  
clear all;
clc;  

%-------基本???置----------------------
thetaT=0;                              %T平台波束斜?角
thetaT=thetaT*pi/180;                  %rad弧度
thetaR=0;                              %R平台波束斜?角
thetaR=thetaR*pi/180;                  %rad弧度
c=2.9979e+8;                           %Speed of light (m/s)光速
fc=1.5e9;                           %Radar center frequency (Hz)?波?率
lambda=c/fc;                           %?波的波?

%%????域
X0=250;         	%方位向[-X0,X0]，范?有自己根据??情?确定
Rtc=3000;       	%?波?射距离
Rrc=3000;       	%?波接收距离
Rc=(Rtc+Rrc)/2; 	%?波距离向
R0=150;         	%距离向[Rc-R0,Rc+R0]

%%距离向??
Tr=1.5e-6;  			%LFM信??? 1.5us (200m)
Br=150e6;     		%LFM信??? 150MHz
Kr=Br/Tr;     		%??斜率
Nr=512;       		%快??采???

%%距离域序列
r=Rc+linspace(-R0,R0,Nr);    %距离域序列
t=2*r/c;                     %距离域t值??
dt=R0*4/c/Nr;                %快??采?周期

%%?率域序列
f=linspace(-1/2/dt,1/2/dt,Nr); %f域序列

%%方位向??
v=6.9621e+03;     %SAR 平台移?速度
Lsar=300;         %合成孔??度
Na=1024;          %慢??采???

%%u域序列
x=linspace(-X0,X0,Na);       %u域序列
u=x/v;                       %u域序列t值??
du=2*X0/v/Na;                %慢??采??隔

%%fu域序列%
fu=linspace(-1/2/du,1/2/du,Na);				%fu域序列
ftdc=v*sin(thetaT);           				%SAR-T平台波束速度
ftdr=-(v*cos(thetaT))^2/lambda/Rtc;
frdc=v*sin(thetaR);           				%SAR-R平台波束速度 
frdr=-(v*cos(thetaR))^2/lambda/Rrc;
fdc=ftdc+frdc;                				%Doppler??中心?率
fdr=ftdr+frdr;                				%Doppler??斜率

%%目?位置
Ntar=3;%目???
Ptar=[  Rrc ,    0 ,  1 		%????：距离向坐?,方位向坐?,sigma             
      Rrc+50,   -50,  1
      Rrc+50,    50,  1];

%%生成SAR正交解?后的回波?据
s_ut=zeros(Nr,Na);  				%生成m×n的double?零矩?，?定?值?存
U=ones(Nr,1)*u;            	%?充?矩?
T=t'*ones(1,Na);           	%?快??拓展?Na列
for i=1:1:Ntar
    rn=Ptar(i,1);          	%目?距离向坐?
    xn=Ptar(i,2);          	%目?方位向坐?
    sigma=Ptar(i,3);       	%目?RCS
    rtn=rn+Rtc-Rrc;        
    RT=sqrt(rtn^2+(rtn*tan(thetaT)+xn-v*U).^2);     %?射目?斜距
    RR=sqrt(rn^2+(rn*tan(thetaT)+xn-v*U).^2);       %接收目?斜距
    R=RT+RR;
    DT=T-R/c;
    phase=pi*Kr*DT.^2-2*pi/lambda*R;
    s_ut=s_ut+sigma*exp(j*phase).*(abs(DT)<Tr/2).*(abs(v*U-xn)<Lsar/2);
end;

%-------距离向??----------------------
%?考信?
p0_t=exp(j*pi*Kr*(t-2*Rc/c).^2).*(abs(t-2*Rc/c)<Tr/2);   %距离向LFM信?
p0_f=fftshift(fft(fftshift(p0_t)));     	% 距离向LFM?考信?的快速傅里???
s_uf=fftshift(fft(fftshift(s_ut))); 			%距离向FFT
src_uf=s_uf.*(conj(p0_f).'*ones(1,Na));  	%距离??
src_ut=fftshift(ifft(fftshift(src_uf))); 	%IFFT后得到距离??后的信?

%-------方位向??----------------------
src_fut=fftshift(fft(fftshift(src_ut).')).'; 		%距离多普勒域

%%二次距离??,距离?移校正原理仿真
src_fuf=fftshift(fft(fftshift(src_uf).')).'; 		%距离??后的二???
F=f'*ones(1,Na);																%?充?矩?
FU=ones(Nr,1)*fu;
p0_2f=exp(j*pi/fc^2/fdr*(FU.*F).^2+j*pi*fdc^2/fc/fdr*F-j*pi/fc/fdr*FU.^2.*F);
s2rc_fuf=src_fuf.*p0_2f;
s2rc_fut=fftshift(ifft(fftshift(s2rc_fuf)));					%距离多普勒域
p0_2fu=exp(j*pi/fdr*(FU-fdc).^2);											%方位向??因子
s2rcac_fut=s2rc_fut.*p0_2fu;													%方位??
s2rcac_fuf=fftshift(fft(fftshift(s2rcac_fut)));				%距离方位??后的二???
s2rcac_ut=fftshift(ifft(fftshift(s2rcac_fut).')).';		%方位向IFFT

%%?果?示
subplot(131)
G=20*log10(abs(s_ut)+1e-6);
gm=max(max(G));
gn=gm-40;																							%?示??范?40dB
G=255/(gm-gn)*(G-gn).*(G>gn);
imagesc(x,r-Rc,-G),colorbar;
grid on,axis tight,
xlabel('Azimuth')
ylabel('Range')
title('Simulated-Signal')

subplot(132)
G=20*log10(abs(src_fut)+1e-6);
gm=max(max(G));
gn=gm-40;																							%?示??范?40dB
G=255/(gm-gn)*(G-gn).*(G>gn);
imagesc(fu,r-Rc,-G),colorbar;
grid on,axis tight,
xlabel('Azimuth')
ylabel('Range')
title('After distance migration correction')

subplot(133)
G=20*log10(abs(s2rc_fut)+1e-6);
gm=max(max(G));
gn=gm-40;																							%?示??范?40dB
G=255/(gm-gn)*(G-gn).*(G>gn);
imagesc(fu,r-Rc,-G),colorbar;
grid on,axis tight,
xlabel('Azimuth')
ylabel('Range')
title('Eliminate phase shift')
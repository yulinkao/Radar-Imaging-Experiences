%�W?�G��_RD��k���p?������u

close all;  
clear all;
clc;  

%-------��???�m----------------------
thetaT=0;                              %T���x�i����?��
thetaT=thetaT*pi/180;                  %rad����
thetaR=0;                              %R���x�i����?��
thetaR=thetaR*pi/180;                  %rad����
c=2.9979e+8;                           %Speed of light (m/s)���t
fc=1.5e9;                           %Radar center frequency (Hz)?�i?�v
lambda=c/fc;                           %?�i���i?

%%????��
X0=250;         	%���V[-X0,X0]�A�S?���ۤv���u??��?�̩w
Rtc=3000;       	%?�i?�g�Z��
Rrc=3000;       	%?�i�����Z��
Rc=(Rtc+Rrc)/2; 	%?�i�Z�æV
R0=150;         	%�Z�æV[Rc-R0,Rc+R0]

%%�Z�æV??
Tr=1.5e-6;  			%LFM�H??? 1.5us (200m)
Br=150e6;     		%LFM�H??? 150MHz
Kr=Br/Tr;     		%??�ײv
Nr=512;       		%��??��???

%%�Z�ð�ǦC
r=Rc+linspace(-R0,R0,Nr);    %�Z�ð�ǦC
t=2*r/c;                     %�Z�ð�t��??
dt=R0*4/c/Nr;                %��??��?�P��

%%?�v��ǦC
f=linspace(-1/2/dt,1/2/dt,Nr); %f��ǦC

%%���V??
v=6.9621e+03;     %SAR ���x��?�t��
Lsar=300;         %�X����??��
Na=1024;          %�C??��???

%%u��ǦC
x=linspace(-X0,X0,Na);       %u��ǦC
u=x/v;                       %u��ǦCt��??
du=2*X0/v/Na;                %�C??��??�j

%%fu��ǦC%
fu=linspace(-1/2/du,1/2/du,Na);				%fu��ǦC
ftdc=v*sin(thetaT);           				%SAR-T���x�i���t��
ftdr=-(v*cos(thetaT))^2/lambda/Rtc;
frdc=v*sin(thetaR);           				%SAR-R���x�i���t�� 
frdr=-(v*cos(thetaR))^2/lambda/Rrc;
fdc=ftdc+frdc;                				%Doppler??����?�v
fdr=ftdr+frdr;                				%Doppler??�ײv

%%��?��m
Ntar=3;%��???
Ptar=[  Rrc ,    0 ,  1 		%????�G�Z�æV��?,���V��?,sigma             
      Rrc+50,   -50,  1
      Rrc+50,    50,  1];

%%�ͦ�SAR�����?�Z���^�i?�u
s_ut=zeros(Nr,Na);  				%�ͦ�m��n��double?�s�x?�A?�w?��?�s
U=ones(Nr,1)*u;            	%?�R?�x?
T=t'*ones(1,Na);           	%?��??�ݮi?Na�C
for i=1:1:Ntar
    rn=Ptar(i,1);          	%��?�Z�æV��?
    xn=Ptar(i,2);          	%��?���V��?
    sigma=Ptar(i,3);       	%��?RCS
    rtn=rn+Rtc-Rrc;        
    RT=sqrt(rtn^2+(rtn*tan(thetaT)+xn-v*U).^2);     %?�g��?�׶Z
    RR=sqrt(rn^2+(rn*tan(thetaT)+xn-v*U).^2);       %������?�׶Z
    R=RT+RR;
    DT=T-R/c;
    phase=pi*Kr*DT.^2-2*pi/lambda*R;
    s_ut=s_ut+sigma*exp(j*phase).*(abs(DT)<Tr/2).*(abs(v*U-xn)<Lsar/2);
end;

%-------�Z�æV??----------------------
%?�ҫH?
p0_t=exp(j*pi*Kr*(t-2*Rc/c).^2).*(abs(t-2*Rc/c)<Tr/2);   %�Z�æVLFM�H?
p0_f=fftshift(fft(fftshift(p0_t)));     	% �Z�æVLFM?�ҫH?���ֳt�Ũ�???
s_uf=fftshift(fft(fftshift(s_ut))); 			%�Z�æVFFT
src_uf=s_uf.*(conj(p0_f).'*ones(1,Na));  	%�Z��??
src_ut=fftshift(ifft(fftshift(src_uf))); 	%IFFT�Z�o��Z��??�Z���H?

%-------���V??----------------------
src_fut=fftshift(fft(fftshift(src_ut).')).'; 		%�Z�æh���ǰ�

%%�G���Z��??,�Z��?���ե���z��u
src_fuf=fftshift(fft(fftshift(src_uf).')).'; 		%�Z��??�Z���G???
F=f'*ones(1,Na);																%?�R?�x?
FU=ones(Nr,1)*fu;
p0_2f=exp(j*pi/fc^2/fdr*(FU.*F).^2+j*pi*fdc^2/fc/fdr*F-j*pi/fc/fdr*FU.^2.*F);
s2rc_fuf=src_fuf.*p0_2f;
s2rc_fut=fftshift(ifft(fftshift(s2rc_fuf)));					%�Z�æh���ǰ�
p0_2fu=exp(j*pi/fdr*(FU-fdc).^2);											%���V??�]�l
s2rcac_fut=s2rc_fut.*p0_2fu;													%���??
s2rcac_fuf=fftshift(fft(fftshift(s2rcac_fut)));				%�Z�ä��??�Z���G???
s2rcac_ut=fftshift(ifft(fftshift(s2rcac_fut).')).';		%���VIFFT

%%?�G?��
subplot(131)
G=20*log10(abs(s_ut)+1e-6);
gm=max(max(G));
gn=gm-40;																							%?��??�S?40dB
G=255/(gm-gn)*(G-gn).*(G>gn);
imagesc(x,r-Rc,-G),colorbar;
grid on,axis tight,
xlabel('Azimuth')
ylabel('Range')
title('Simulated-Signal')

subplot(132)
G=20*log10(abs(src_fut)+1e-6);
gm=max(max(G));
gn=gm-40;																							%?��??�S?40dB
G=255/(gm-gn)*(G-gn).*(G>gn);
imagesc(fu,r-Rc,-G),colorbar;
grid on,axis tight,
xlabel('Azimuth')
ylabel('Range')
title('After distance migration correction')

subplot(133)
G=20*log10(abs(s2rc_fut)+1e-6);
gm=max(max(G));
gn=gm-40;																							%?��??�S?40dB
G=255/(gm-gn)*(G-gn).*(G>gn);
imagesc(fu,r-Rc,-G),colorbar;
grid on,axis tight,
xlabel('Azimuth')
ylabel('Range')
title('Eliminate phase shift')
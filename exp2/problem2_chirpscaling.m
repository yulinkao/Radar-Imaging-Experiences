%名?：基于chirp_scaling算法的雷?成像仿真

close all;  
clear all;  
clc;

%-------基本???置----------------------
Tr=200;							%??200m
Br=1;								%??
Kr=Br/Tr;						%??斜率
Fc=4;								%??
Nr=512;							%快??采???
Na=1024;						%慢??采???

%%????域
Xc=1200;X0=250;										%定?距离向范?
x=Xc+linspace(-X0,X0,Nr);					%x域序列:Xc-X0~Xc+X0
dx=2*X0/Nr;												%定?步?
kx=linspace(-1/dx/2,1/dx/2,Nr);		%kx域序列
Y0=200;
y=linspace(-Y0,Y0,Na);						%y域序列:-Y0~Y0
dy=2*Y0/Na;
ky=linspace(-1/dy/2,1/dy/2,Na);		%ky域序列

%%方位向??
Lsar=300;				%??300m,合成孔??度
Ba=1;						%??1(1/m)
Ka=Fc/Xc;				%??斜率 Ka=Ba/Ta=Fc/Xc

%%目?几何?系target geometry
%x坐?,y坐?,复后向散射系? 
Ptar=[Xc,0,1+0j              
      Xc+50,-50,1+0j
      Xc+50,50,1+0j
      Xc-50,-50,1+0j
      Xc-50,50,1+0j];

%%生成SAR正交解?后的回波?据
Srnm=zeros(Nr,Na);
N=size(Ptar,1);																		%目???
for i=1:1:N
    xn=Ptar(i,1);yn=Ptar(i,2);sigma=Ptar(i,3);		%提取每?目?的信息
    X=x.'*ones(1,Na);															%?充?矩?
    Y=ones(Nr,1)*y;																%?充?矩?
    DX=sqrt(xn^2+(Y-yn).^2);											%中??量
    phase=pi*Kr*(X-DX).^2-2*pi*Fc*DX;							%回波相位
    
    %回波累加
    Srnm=Srnm+sigma*exp(j*phase).*(abs(X-DX)<Tr/2).*(abs(Y-yn)<Lsar/2);
end

%%?据准?
phi0=-x'*sqrt(Fc^2-ky.^2);
phi1=-Fc*x'*(1./sqrt(Fc^2-ky.^2));
phi2=1/2*x'*(ky.^2./(Fc^2-ky.^2).^1.5);
Cs=ones(Nr,1)*(Fc./sqrt(Fc^2-ky.^2)-1);
Ks=1./(1/Kr-2*phi2);

%-------方位向??----------------------
s_xky=fftshift(fft(fftshift(Srnm).')).';						%方位向FFT

%-------Chirp Scaling-----------------
scs_xky=s_xky.*exp(j*pi*Cs.*Ks.*(x'*ones(1,Na)-Xc*(1+Cs)).^2);
s1=ifft(scs_xky);%??示存??据

%-------距离向??----------------------
scs_kxky=fftshift(fft(fftshift(scs_xky)));					%距离向FFT

%%距离?移校正&距离向匹配?波
srmc_kxky=scs_kxky.*exp(j*pi*(kx.^2'*ones(1,Na))./(1+Cs)./Ks...
                     +j*2*pi*Xc*Cs.*(kx'*ones(1,Na)));
srmc_xky=fftshift(ifft(fftshift(srmc_kxky)));				%距离向IFFT

%%消除?差函?.方位向匹配?波
f_xky=srmc_xky.*exp(-j*pi*Ks.*Cs.*(1+Cs).*((x-Xc).^2'*ones(1,Na))...
           -j*2*pi*phi0);
f_xy=fftshift(ifft(fftshift(f_xky).')).';						%方位向IFFT

%%?s1?行距离向???示
p0_x=exp(j*pi*Kr*(x-Xc).^2).*(abs(x-Xc)<Tr/2);			%距离向LFM信?
p0_kx=fftshift(fft(fftshift(p0_x)));
p0_y=exp(-j*pi*Ka*y.^2).*(abs(y)<Lsar/2);						%方位向LFM信?
p0_ky=fftshift(fft(fftshift(p0_y))); 
s_kxy=fftshift(fft(fftshift(s1)));									%距离向FFT
sxc_kxy=s_kxy.*(conj(p0_kx).'*ones(1,Na));
sxc_kxky=fftshift(fft(fftshift(sxc_kxy).')).';			%距离??后的2D?域信?
sxc_xy=fftshift(ifft(fftshift(sxc_kxy)));						%距离??后的信?
sxc_xky=fftshift(fft(fftshift(sxc_xy).')).';				%距离??后.距离-多普勒域

%%?果?示
subplot(131)
colorbar;
imagesc(255-abs(Srnm));
grid on,axis tight,
xlabel('Azimuth'),ylabel('Range'),
title('Simulated-Signal');

subplot(132)
colorbar;
imagesc(255-abs(srmc_xky));
grid on,axis tight,
xlabel('Azimuth'),ylabel('Range'),
title('After distance migration correction');

subplot(133)
colorbar;
imagesc(255-abs(f_xky)); 
grid on,axis tight,
xlabel('Azimuth'),ylabel('Range'),
title('Eliminate phase shift');

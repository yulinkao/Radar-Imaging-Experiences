%�W?�G��_chirp_scaling��k���p?������u

close all;  
clear all;  
clc;

%-------��???�m----------------------
Tr=200;							%??200m
Br=1;								%??
Kr=Br/Tr;						%??�ײv
Fc=4;								%??
Nr=512;							%��??��???
Na=1024;						%�C??��???

%%????��
Xc=1200;X0=250;										%�w?�Z�æV�S?
x=Xc+linspace(-X0,X0,Nr);					%x��ǦC:Xc-X0~Xc+X0
dx=2*X0/Nr;												%�w?�B?
kx=linspace(-1/dx/2,1/dx/2,Nr);		%kx��ǦC
Y0=200;
y=linspace(-Y0,Y0,Na);						%y��ǦC:-Y0~Y0
dy=2*Y0/Na;
ky=linspace(-1/dy/2,1/dy/2,Na);		%ky��ǦC

%%���V??
Lsar=300;				%??300m,�X����??��
Ba=1;						%??1(1/m)
Ka=Fc/Xc;				%??�ײv Ka=Ba/Ta=Fc/Xc

%%��?�L��?�ttarget geometry
%x��?,y��?,�`�Z�V���g�t? 
Ptar=[Xc,0,1+0j              
      Xc+50,-50,1+0j
      Xc+50,50,1+0j
      Xc-50,-50,1+0j
      Xc-50,50,1+0j];

%%�ͦ�SAR�����?�Z���^�i?�u
Srnm=zeros(Nr,Na);
N=size(Ptar,1);																		%��???
for i=1:1:N
    xn=Ptar(i,1);yn=Ptar(i,2);sigma=Ptar(i,3);		%�����C?��?���H��
    X=x.'*ones(1,Na);															%?�R?�x?
    Y=ones(Nr,1)*y;																%?�R?�x?
    DX=sqrt(xn^2+(Y-yn).^2);											%��??�q
    phase=pi*Kr*(X-DX).^2-2*pi*Fc*DX;							%�^�i�ۦ�
    
    %�^�i�֥[
    Srnm=Srnm+sigma*exp(j*phase).*(abs(X-DX)<Tr/2).*(abs(Y-yn)<Lsar/2);
end

%%?�u��?
phi0=-x'*sqrt(Fc^2-ky.^2);
phi1=-Fc*x'*(1./sqrt(Fc^2-ky.^2));
phi2=1/2*x'*(ky.^2./(Fc^2-ky.^2).^1.5);
Cs=ones(Nr,1)*(Fc./sqrt(Fc^2-ky.^2)-1);
Ks=1./(1/Kr-2*phi2);

%-------���V??----------------------
s_xky=fftshift(fft(fftshift(Srnm).')).';						%���VFFT

%-------Chirp Scaling-----------------
scs_xky=s_xky.*exp(j*pi*Cs.*Ks.*(x'*ones(1,Na)-Xc*(1+Cs)).^2);
s1=ifft(scs_xky);%??�ܦs??�u

%-------�Z�æV??----------------------
scs_kxky=fftshift(fft(fftshift(scs_xky)));					%�Z�æVFFT

%%�Z��?���ե�&�Z�æV�ǰt?�i
srmc_kxky=scs_kxky.*exp(j*pi*(kx.^2'*ones(1,Na))./(1+Cs)./Ks...
                     +j*2*pi*Xc*Cs.*(kx'*ones(1,Na)));
srmc_xky=fftshift(ifft(fftshift(srmc_kxky)));				%�Z�æVIFFT

%%����?�t��?.���V�ǰt?�i
f_xky=srmc_xky.*exp(-j*pi*Ks.*Cs.*(1+Cs).*((x-Xc).^2'*ones(1,Na))...
           -j*2*pi*phi0);
f_xy=fftshift(ifft(fftshift(f_xky).')).';						%���VIFFT

%%?s1?��Z�æV???��
p0_x=exp(j*pi*Kr*(x-Xc).^2).*(abs(x-Xc)<Tr/2);			%�Z�æVLFM�H?
p0_kx=fftshift(fft(fftshift(p0_x)));
p0_y=exp(-j*pi*Ka*y.^2).*(abs(y)<Lsar/2);						%���VLFM�H?
p0_ky=fftshift(fft(fftshift(p0_y))); 
s_kxy=fftshift(fft(fftshift(s1)));									%�Z�æVFFT
sxc_kxy=s_kxy.*(conj(p0_kx).'*ones(1,Na));
sxc_kxky=fftshift(fft(fftshift(sxc_kxy).')).';			%�Z��??�Z��2D?��H?
sxc_xy=fftshift(ifft(fftshift(sxc_kxy)));						%�Z��??�Z���H?
sxc_xky=fftshift(fft(fftshift(sxc_xy).')).';				%�Z��??�Z.�Z��-�h���ǰ�

%%?�G?��
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

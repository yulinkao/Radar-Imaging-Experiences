%�W?�G��_???�u?�z���p?������u
%% ???�u?�z�A��?��Z��??�A�M�Z��Z�è�???�A�̦Z��즨��

close all;
clear all;
clc;

%-------?�u?��-------------------------
load data_Raw.mat; 
[Nrn,Nan]=size(data1);  			%�Z��????�M���????

%-------��???�m----------------------
%% ?�p?�t???�t�m
%%�`?�w?
C=2.9979e+8;                  %���t
%%�p???
Fc=5.300e+9;                  %??
lambda=C/Fc;                  %�i?

%%��??��??
Xmin=-150;                    %��??����V�S?[Xmin,Xmax]
Xmax=150;                  		%����?�줤?
Yc=8000;                      %��??��Z�æV�S?[Yc-Y0,Yc+Y0]	
Y0=1000;                      %����?��?2*Y0
                              
%%?�D??
V=6.9621e+03;                 %SAR��??�t��
H=6000;                       %����
R0=sqrt(Yc^2+H^2);            %�̵u�Z��

%%��???
D=4;                          %���V��??��
Lsar=lambda*R0/D;   					%SAR�X����??��
Tsar=Lsar/V;                  %SAR�Ӯg??

%%�C??��??
Ka=-2*V^2/lambda/R0;    						%�h����?��??�v
Ba=abs(Ka*Tsar);           					%�h����?�v?��??
PRF=1256.98;            						%??���`?�v.�S?�`?�v
PRT=1/PRF;                  				%??���`??
ds=PRT;                       			%�C?�쪺??�B?
Nslow=ceil((Xmax-Xmin+Lsar)/V/ds); 	%�C?�쪺��??
Nslow=2^nextpow2(Nslow);     				%nextpow2?�̾a��2��?����?.?fft??����?
sn=linspace((Xmin-Lsar/2)/V,(Xmax+Lsar/2)/V,Nslow);			%�C??�쪺??�x?
PRT=(Xmax-Xmin+Lsar)/V/Nslow;
PRF=1/PRT;
ds=PRT;

%%��??��???�m
Tr=0.5e-6;                    %??��???5us
Br=100e6;                     %chirp?�v?��???30MHz
Kr=Br/Tr;                     %chirp??�v
Fsr=2*Br;                     %��?���??�v.?3����??
dt=1/Fsr;                     %��?���??�j
Rmin=sqrt((Yc-Y0)^2+H^2);
Rmax=sqrt((Yc+Y0)^2+H^2+(Lsar/2)^2);
Nfast=ceil(2*(Rmax-Rmin)/C/dt+Tr/dt);				%��?�쪺��??�q
Nfast=2^nextpow2(Nfast);                   	%��s?2��?��.��K?��fft??
tm=linspace(2*Rmin/C,2*Rmax/C+Tr,Nfast); 		%��?�쪺�ô�??�x?
dt=(2*Rmax/C+Tr-2*Rmin/C)/Nfast;    				%��s?�j
Fsr=1/dt;

%%����v???�m
DY=C/2/Br;                    %�Z�æV����v
DX=D/2;                       %���V����v

%%�ͦ��^�i�H?
K=Nrn;                                		%��??��
N=Nslow;                                  %�C?�쪺��??
M=Nfast;                                  %��?�쪺��??
T=data1;                                	%��?�x?
Srnm=zeros(N,M);                          %�ͦ��s�x?�s?�^�i�H?
for k=1:1:K                               %?�@K?��?
    sigma=T(k,3);                         %�o���?���Ϯg�t?
    Dslow=sn*V-T(k,1);                    %���V�Z��.��v����V���Z��
    R=sqrt(Dslow.^2+T(k,2)^2+H^2);        %??�Z�ïx?
    tau=2*R/C;                            %�^�i��?�_?�g�i����?
    Dfast=ones(N,1)*tm-tau'*ones(1,M);    
    phase=pi*Kr*Dfast.^2-(4*pi/lambda)*(R'*ones(1,M));		%�ۦ�
    Srnm=Srnm+sigma*exp(j*phase).*(0<Dfast&Dfast<Tr).*((abs(Dslow)<Lsar/2)'*ones(1,M));
end

%-------�Z�æV??----------------------
tic;
tr=tm-2*Rmin/C;
Refr=exp(j*pi*Kr*tr.^2).*(0<tr&tr<Tr);
Sr=ifft(fft(Srnm).*(ones(N,1)*conj(fft(Refr))));
Gr=abs(Sr);
Sa_RD = fft(Sr);     									%���VFFT??�Z�æh����?��Z��?���ե�

Kp=1;                        					%?��Ϊ�????�i��

%%����?��Z��?���q�x?
for n=1:N     												%?�@��N?����?
    for m=1:M 												%�C?����?�W��M?�Z�ê�?
        delta_R = (1/8)*(lambda/V)^2*(R0+(m-M/2)*C/2/Fsr)*((n-N/2)*PRF/N)^2;
        RMC=2*delta_R*Fsr/C;    			%��??delta_R/DY.�Z�î{?�F�L?�Z��?��
        delta_RMC = RMC-round(RMC);		%�Z�î{?�q���p?����
        if m+round(RMC)>M             %�P?�O�_�W�X?��
            Sa_RD(n,m)=Sa_RD(n,M/2);   
        else
            if delta_RMC>=0.5  				%���J
                Sa_RD(n,m)=Sa_RD(n,m+round(RMC)+1);
            else               				%�|��
                Sa_RD(n,m)=Sa_RD(n,m+round(RMC));
            end
        end
    end
end

%-------�Z�ñp?�ե�--------------------
Sr_rmc=ifft(Sa_RD);   
Ga = abs(Sr_rmc);

%-------���V??----------------------
ta=sn-Xmin/V;
Refa=exp(j*pi*Ka*ta.^2).*(abs(ta)<Tsar/2);
Sa=ifft(fft(Sr_rmc).*(conj(fft(Refa)).'*ones(1,M)));
Gar=abs(Sa);
toc;

%%?�G?��
figure(1)
subplot(121);
colorbar;
row=tm*C/2-2008;col=sn*V-26;
imagesc(row,col,255-Gr);           						%�Z�æV??.���ե��Z�ñp?��?��
axis([Yc-Y0,Yc+Y0,Xmin-Lsar/2,Xmax+Lsar/2]);
xlabel('Azimuth');ylabel('Range');
title('Azumuth compress without RCMC');

subplot(122);
colorbar;
imagesc(row,col,255-Ga);          						%�Z�æV??.�ե��Z�ñp?�Z��?��
axis([Yc-Y0,Yc+Y0,Xmin-Lsar/2,Xmax+Lsar/2]);
xlabel('Azimuth');ylabel('Range');
title('Azumuth compress with RCMC');

figure(2)
colorbar;
imagesc(row,col,255-Gar);          						%���V??�Z��?��
axis([Yc-Y0,Yc+Y0,Xmin-Lsar/2,Xmax+Lsar/2]);
xlabel('Azimuth');ylabel('Range');
title('After Range Compress');

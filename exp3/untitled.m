%% ???�u?�z�A��?��Z��??�A�M�Z��Z�è�???�A�̦Z��즨��
clear all;clc;close all;
%% ?�u?��
load data_Raw.mat; 
[Nrn,Nan]=size(data1);  % �Z��????�M���????
%%�`?�w?
C=2.9979e+8;                            %���t
%%�p???
Fc=1.5e9;                          %??1GHz
lambda=C/Fc;                     %�i?
%%��??��??
Xmin=-200;                          %��??����V�S?[Xmin,Xmax]
Xmax=200;                  
Yc=8000;                      %����?�줤?
Y0=1000;                          %��??��Z�æV�S?[Yc-Y0,Yc+Y0]
                                       %����?��?2*Y0
%%?�D??
V=120;                            %SAR��??�t��120 m/s
H=6000;                          %���� 6000 m
R0=sqrt(Yc^2+H^2);               %�̵u�Z��
%%��???
D=4;                                %���V��??��
Lsar=lambda*R0/D;   %SAR�X����??��.�m�X����?�p?�����X�X��k�O??�nP.100
Tsar=Lsar/V;                   %SAR�Ӯg??
%%�C??��??
Ka=-2*V^2/lambda/R0;    %�h����?��??�vP.93
Ba=abs(Ka*Tsar);           %�h����?�v?��??
PRF=Ba;                         %??���`?�v.PRF��??�h����?�v����?�v.�S?�`?�v.�ҥH���_Ba.P.93
PRT=1/PRF;                   %??���`??
ds=PRT;                         %�C?�쪺??�B?
Nslow=ceil((Xmax-Xmin+Lsar)/V/ds); %�C?�쪺��??.ceil?�����?.?�XP.76��?�z��
Nslow=2^nextpow2(Nslow);              %nextpow2?�̾a��2��?����?.?��?fft??����?
sn=linspace((Xmin-Lsar/2)/V,(Xmax+Lsar/2)/V,Nslow);%�C??�쪺??�x?
PRT=(Xmax-Xmin+Lsar)/V/Nslow;    %�Ѥ_Nslow��?�F.�ҥH��?���@��??�]�ݭn��s.�P��?�p�F
PRF=1/PRT;
ds=PRT;
%%��??��???�m
Tr=0.5e-6;                         %??��???5us
Br=100e6;                        %chirp?�v?��???30MHz
Kr=Br/Tr;                        %chirp??�v
Fsr=2*Br;                        %��?���??�v.?3����??
dt=1/Fsr;                         %��?���??�j
Rmin=sqrt((Yc-Y0)^2+H^2);
Rmax=sqrt((Yc+Y0)^2+H^2+(Lsar/2)^2);
Nfast=ceil(2*(Rmax-Rmin)/C/dt+Tr/dt);%��?�쪺��??�q
Nfast=2^nextpow2(Nfast);                   %��s?2��?��.��K?��fft??
tm=linspace(2*Rmin/C,2*Rmax/C+Tr,Nfast); %��?�쪺�ô�??�x?
dt=(2*Rmax/C+Tr-2*Rmin/C)/Nfast;    %��s?�j
Fsr=1/dt;
%%����v???�m
DY=C/2/Br;                           %�Z�æV����v
DX=D/2;                                %���V����v


%%================================================================
%%�ͦ��^�i�H?
K=Nrn;                                %��??��
N=Nslow;                                  %�C?�쪺��??
M=Nfast;                                  %��?�쪺��??
T=data1;                                %��?�x?
Srnm=zeros(N,M);                          %�ͦ��s�x?�s?�^�i�H?
for k=1:1:K                               %?�@K?��?
    sigma=T(k,3);                         %�o���?���Ϯg�t?
    Dslow=sn*V-T(k,1);                    %���V�Z��.��v����V���Z��
    R=sqrt(Dslow.^2+T(k,2)^2+H^2);        %??�Z�ïx?
    tau=2*R/C;                            %�^�i��?�_?�g�i����?
    Dfast=ones(N,1)*tm-tau'*ones(1,M);    %(t-tau).��?�N�O??�x?.ones(N,1)�Mones(1,M)���O?�F?��?�i?�x?
    phase=pi*Kr*Dfast.^2-(4*pi/lambda)*(R'*ones(1,M));%�ۦ�.����??P.96
    Srnm=Srnm+sigma*exp(j*phase).*(0<Dfast&Dfast<Tr).*((abs(Dslow)<Lsar/2)'*ones(1,M));%�Ѥ_�O�h?��?�Ϯg���^�i.�ҥH��??��?�[
end
%%================================================================
%%�Z��-�h���Ǻ�k?�l
%%�Z�æV??

tr=tm-2*Rmin/C;
Refr=exp(j*pi*Kr*tr.^2).*(0<tr&tr<Tr);
Sr=ifft(fft(Srnm).*(ones(N,1)*conj(fft(Refr))));
Gr=abs(Sr);
%%?�l?��Z��?��??��???���Z�è�?? �D�n�O�]?�׶Z��?�Ƥް_�^�i�]?���p? 
%%??��k�G�̪�?�촡�Ȫk.���^?�G��??��Z�æh���ǰ�.��????����??��X�Z�ñp?�q.�o��Z�ñp?�q�O�Z�ä���v�����.
%%?��ȥi��?�p?.���ӥ|�٤��J����k���?��?.�ӦZ�b?����?�W?�h�p?�q
%%���V��fft?�z �A�b?�찵�Z��?��??
Sa_RD = fft(Sr);     %  ���VFFT ??�Z�æh����?��Z��?���ե�
%�Z�ñp??��,�Ѥ_�O��?? .fdc=0,�u�ݭn?��Z��?��??�C
Kp=1;                                  %?��Ϊ�????�i��

%%����?��Z��?���q�x?
for n=1:N     %?�@��N?����?
    for m=1:M %�C?����?�W��M?�Z�ê�?
        delta_R = (1/8)*(lambda/V)^2*(R0+(m-M/2)*C/2/Fsr)*((n-N/2)*PRF/N)^2;%�Z��?���qP.160�F(R0+(m-M/2)*C/2/Fsr)�G�C?�Z�æV?m��R0��s�F(n-N/2)*PRF/N�G���P���V���h����?�v���@?
        RMC=2*delta_R*Fsr/C;    %��??delta_R/DY.�Z�î{?�F�L?�Z��?��
        delta_RMC = RMC-round(RMC);%�Z�î{?�q���p?����
        if m+round(RMC)>M              %�P?�O�_�W�X?��
            Sa_RD(n,m)=Sa_RD(n,M/2);   
        else
            if delta_RMC>=0.5  %���J
                Sa_RD(n,m)=Sa_RD(n,m+round(RMC)+1);
            else               %�|��
                Sa_RD(n,m)=Sa_RD(n,m+round(RMC));
            end
        end
    end
end

%========================
Sr_rmc=ifft(Sa_RD);   %%�Z�ñp?�ե��Z?���?��
Ga = abs(Sr_rmc);
%%���V??
ta=sn-Xmin/V;
Refa=exp(j*pi*Ka*ta.^2).*(abs(ta)<Tsar/2);
Sa=ifft(fft(Sr_rmc).*(conj(fft(Refa)).'*ones(1,M)));
Gar=abs(Sa);

%%================================================================
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


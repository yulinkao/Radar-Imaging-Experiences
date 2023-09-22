%% ???据?理，先?行距离??，然后把距离走???，最后方位成像
clear all;clc;close all;
%% ?据?取
load data_Raw.mat; 
[Nrn,Nan]=size(data1);  % 距离????和方位????
%%常?定?
C=2.9979e+8;                            %光速
%%雷???
Fc=1.5e9;                          %??1GHz
lambda=C/Fc;                     %波?
%%目??域??
Xmin=-200;                          %目??域方位向范?[Xmin,Xmax]
Xmax=200;                  
Yc=8000;                      %成像?域中?
Y0=1000;                          %目??域距离向范?[Yc-Y0,Yc+Y0]
                                       %成像?度?2*Y0
%%?道??
V=120;                            %SAR的??速度120 m/s
H=6000;                          %高度 6000 m
R0=sqrt(Yc^2+H^2);               %最短距离
%%天???
D=4;                                %方位向天??度
Lsar=lambda*R0/D;   %SAR合成孔??度.《合成孔?雷?成像——算法与??》P.100
Tsar=Lsar/V;                   %SAR照射??
%%慢??域??
Ka=-2*V^2/lambda/R0;    %多普勒?域??率P.93
Ba=abs(Ka*Tsar);           %多普勒?率?制??
PRF=Ba;                         %??重复?率.PRF其??多普勒?率的采?率.又?复?率.所以等于Ba.P.93
PRT=1/PRF;                   %??重复??
ds=PRT;                         %慢?域的??步?
Nslow=ceil((Xmax-Xmin+Lsar)/V/ds); %慢?域的采??.ceil?取整函?.?合P.76的?理解
Nslow=2^nextpow2(Nslow);              %nextpow2?最靠近2的?次函?.?里?fft??做准?
sn=linspace((Xmin-Lsar/2)/V,(Xmax+Lsar/2)/V,Nslow);%慢??域的??矩?
PRT=(Xmax-Xmin+Lsar)/V/Nslow;    %由于Nslow改?了.所以相?的一些??也需要更新.周期?小了
PRF=1/PRT;
ds=PRT;
%%快??域???置
Tr=0.5e-6;                         %??持???5us
Br=100e6;                        %chirp?率?制???30MHz
Kr=Br/Tr;                        %chirp??率
Fsr=2*Br;                        %快?域采??率.?3倍的??
dt=1/Fsr;                         %快?域采??隔
Rmin=sqrt((Yc-Y0)^2+H^2);
Rmax=sqrt((Yc+Y0)^2+H^2+(Lsar/2)^2);
Nfast=ceil(2*(Rmax-Rmin)/C/dt+Tr/dt);%快?域的采??量
Nfast=2^nextpow2(Nfast);                   %更新?2的?次.方便?行fft??
tm=linspace(2*Rmin/C,2*Rmax/C+Tr,Nfast); %快?域的离散??矩?
dt=(2*Rmax/C+Tr-2*Rmin/C)/Nfast;    %更新?隔
Fsr=1/dt;
%%分辨率???置
DY=C/2/Br;                           %距离向分辨率
DX=D/2;                                %方位向分辨率


%%================================================================
%%生成回波信?
K=Nrn;                                %目??目
N=Nslow;                                  %慢?域的采??
M=Nfast;                                  %快?域的采??
T=data1;                                %目?矩?
Srnm=zeros(N,M);                          %生成零矩?存?回波信?
for k=1:1:K                               %?共K?目?
    sigma=T(k,3);                         %得到目?的反射系?
    Dslow=sn*V-T(k,1);                    %方位向距离.投影到方位向的距离
    R=sqrt(Dslow.^2+T(k,2)^2+H^2);        %??距离矩?
    tau=2*R/C;                            %回波相?于?射波的延?
    Dfast=ones(N,1)*tm-tau'*ones(1,M);    %(t-tau).其?就是??矩?.ones(N,1)和ones(1,M)都是?了?其?展?矩?
    phase=pi*Kr*Dfast.^2-(4*pi/lambda)*(R'*ones(1,M));%相位.公式??P.96
    Srnm=Srnm+sigma*exp(j*phase).*(0<Dfast&Dfast<Tr).*((abs(Dslow)<Lsar/2)'*ones(1,M));%由于是多?目?反射的回波.所以此??行?加
end
%%================================================================
%%距离-多普勒算法?始
%%距离向??

tr=tm-2*Rmin/C;
Refr=exp(j*pi*Kr*tr.^2).*(0<tr&tr<Tr);
Sr=ifft(fft(Srnm).*(ones(N,1)*conj(fft(Refr))));
Gr=abs(Sr);
%%?始?行距离?曲??正???有距离走?? 主要是因?斜距的?化引起回波包?的徙? 
%%??方法：最近?域插值法.具体?：先??到距离多普勒域.分????像素??算出距离徙?量.得到距离徙?量与距离分辨率的比值.
%%?比值可能?小?.按照四舍五入的方法近似?整?.而后在?像素?上?去徙?量
%%方位向做fft?理 再在?域做距离?曲??
Sa_RD = fft(Sr);     %  方位向FFT ??距离多普域?行距离?曲校正
%距离徙??算,由于是正?? .fdc=0,只需要?行距离?曲??。
Kp=1;                                  %?算或者????波比

%%首先?算距离?移量矩?
for n=1:N     %?共有N?方位采?
    for m=1:M %每?方位采?上有M?距离采?
        delta_R = (1/8)*(lambda/V)^2*(R0+(m-M/2)*C/2/Fsr)*((n-N/2)*PRF/N)^2;%距离?移量P.160；(R0+(m-M/2)*C/2/Fsr)：每?距离向?m的R0更新；(n-N/2)*PRF/N：不同方位向的多普勒?率不一?
        RMC=2*delta_R*Fsr/C;    %此??delta_R/DY.距离徒?了几?距离?元
        delta_RMC = RMC-round(RMC);%距离徒?量的小?部分
        if m+round(RMC)>M              %判?是否超出?界
            Sa_RD(n,m)=Sa_RD(n,M/2);   
        else
            if delta_RMC>=0.5  %五入
                Sa_RD(n,m)=Sa_RD(n,m+round(RMC)+1);
            else               %四舍
                Sa_RD(n,m)=Sa_RD(n,m+round(RMC));
            end
        end
    end
end

%========================
Sr_rmc=ifft(Sa_RD);   %%距离徙?校正后?原到?域
Ga = abs(Sr_rmc);
%%方位向??
ta=sn-Xmin/V;
Refa=exp(j*pi*Ka*ta.^2).*(abs(ta)<Tsar/2);
Sa=ifft(fft(Sr_rmc).*(conj(fft(Refa)).'*ones(1,M)));
Gar=abs(Sa);

%%================================================================
%%?果?示
figure(1)
subplot(121);
colorbar;
row=tm*C/2-2008;col=sn*V-26;
imagesc(row,col,255-Gr);           						%距离向??.未校正距离徙?的?像
axis([Yc-Y0,Yc+Y0,Xmin-Lsar/2,Xmax+Lsar/2]);
xlabel('Azimuth');ylabel('Range');
title('Azumuth compress without RCMC');

subplot(122);
colorbar;
imagesc(row,col,255-Ga);          						%距离向??.校正距离徙?后的?像
axis([Yc-Y0,Yc+Y0,Xmin-Lsar/2,Xmax+Lsar/2]);
xlabel('Azimuth');ylabel('Range');
title('Azumuth compress with RCMC');

figure(2)
colorbar;
imagesc(row,col,255-Gar);          						%方位向??后的?像
axis([Yc-Y0,Yc+Y0,Xmin-Lsar/2,Xmax+Lsar/2]);
xlabel('Azimuth');ylabel('Range');
title('After Range Compress');


clear all
close all
clc 

%filterDesigner
%%
%filtr=designfilt
%%
fvtool(filtr)
%% FIR
fvz = 100;
a=1;
n=64;
%b=fir1(n,20/(fvz/2),'low')
%b=fir1(n,20/(fvz/2),'high')
b=fir1(n,[20/(fvz/2) 30/(fvz/2)],'bandpass')
zplane(b,a)
figure 
freqz(b,a)
figure 
plot(abs(fft(b,fvz)))

%% IIR
fvz = 100;
n=16;
[b,a] = butter(n,20/(fvz/2),'low')
zplane(b,a)
figure 
freqz(b,a)
figure 
plot(abs(fft(impz(b,a),fvz)))

%%
%filterDesigner

%%
load('signal_3.mat')
load('DP_FIR.mat')
b = DP_FIR.Numerator;
a=1;
y1=filter(b,a,x)

y2 = conv(x,b,'full')

figure 
subplot 311
plot(x)
xlim([0,2800])
subplot 312
plot(y1)
xlim([0,2800])
subplot 313
plot(y2)
xlim([0,2800])

%%
load('signal_3.mat')
load('DP_FIR.mat')
b = DP_FIR.Numerator;
a=1;

%doplnění nulami
x1=[zeros(1,length(b)), x, zeros(1,length(b))];

%doplnění první a poslední hodnotou
x2=[ones(1,length(b))*x(1), x, ones(1,length(b))*x(end)];
plot(x2)


y1 = conv(x1,b,'full')
y1 = y1(length(b):end-length(b))

y2 = conv(x1,b,'full')
y2 = y2(length(b):end-length(b))

figure 
subplot 311
plot(x)
xlim([0,2800])
subplot 312
plot(y1)
xlim([0,2800])
subplot 313
plot(y2)
xlim([0,2800])

%% ideální filtr 
filterDesigner

%%
load('signal_3.mat')
load('DP_IIR.mat')
SOS=DP_IIR.SOSMatrix;
y1=sosfilt(SOS,x)
plot(y1)

%%
y2=filter(DP_IIR,x)
figure
subplot 311
plot(x)
xlim([0,2800])
subplot 312
plot(y1)
xlim([0,2800])
subplot 313
plot(y2)
xlim([0,2800])
%% 
clear all
close all
clc


%% FIR pasmova zadrz 
load('signal_3.mat')
load('FIR_PZ.mat')
load('IIR_PZ_II.mat')

lin=filter(FIR_PZ, x);
nelin=filter(IIR_PZ_II, x);
plot(lin)

%%
lin=lin((length(FIR_PZ.Numerator)-1)/2:end);

%%

figure
subplot 311
plot(x)
xlim([0,1200])
subplot 312
plot(lin)
xlim([0,1200])
subplot 313
plot(nelin)
xlim([0,1200])

%%
load('signal_orig.mat')
x = [0,0,0,x];

figure 
plot(x,'g')
hold on 
plot (lin,'b')
plot(nelin, 'r')
xlim([1950,2100])
legend('orig','FIR', 'IIR')

%% 

b_rek = [1,0,0,0,0,-1]
a_rek = [1,-1]
figure 
zplane(b_rek, a_rek)

b_nerek = [1,1,1,1,1]
a_nerek = [1]
figure 
zplane(b_nerek, a_nerek)
%%
clear all 
close all
clc

%% mediánové filtry 
%x = [0,0,0,1,0,0,0];
%x= [0,0,0,1,1,0,0,0];
x=[0,0,0,1,1,1,0,0,0];

delka_okna=7
zpozdeni=(delka_okna-1)/2
y=[]
y2=[]

for i=1:length(x)-delka_okna+1
    y(i)=median(x(i:i+delka_okna-1))
    y2(i)=mean(x(i:i+delka_okna-1))
end 

figure 
subplot 311
stem(x)
subplot 312
stem([zeros(1,zpozdeni) y zeros(1,zpozdeni)])
subplot 313
stem([zeros(1,zpozdeni) y2 zeros(1,zpozdeni)])

%% prumerujici filtr 
clear all
close all 
clc

load ('ekg_250.mat')
x_sum= impulse_noise(x,0.2,1)


delka_okna=15
zpozdeni=(delka_okna-1)/2
y=[]
y2=[]

for i=1:length(x_sum)-delka_okna+1
    y(i)=median(x_sum(i:i+delka_okna-1))
    y2(i)=mean(x_sum(i:i+delka_okna-1))
end 

figure 
subplot 221
plot(x)
subplot 222
plot(x_sum)
subplot 223
plot([zeros(1,zpozdeni) y zeros(1,zpozdeni)])
title('medfilt')
subplot 224
plot([zeros(1,zpozdeni) y2 zeros(1,zpozdeni)])
title('mean')

%% kumulace 
clear all 
close all
clc 

x = [8,5,6,1,8,7,5,7,2,7,7,2,6,7,3,6,8,8,2,7,7];
perioda =3;
delka_okna =3;
Mmax=length(x)/perioda;
R=[];
for i=1:Mmax
    R(:,i)=x((i-1)*perioda+1:(i-1)*perioda+perioda)
end 

%% kumulace s pevnym oknem 
kum=zeros(perioda,1)
for i=1:Mmax
    kum=1/Mmax*R(:,i)+kum;
    figure
    stem(kum)
    xlim([0,4])
    ylim([0,10])
    pause(1)
end 

%% kumulace s plovoucim oknem 

perioda =3;
delka_okna =3;
Mmax=length(x)/perioda;
kum=zeros(perioda,1)

for i=1:Mmax
    if(i-delka_okna)<=0
        kum=1/delka_okna*R(:,i)+kum;
    else 
        kum=1/delka_okna*R(:,i)+kum-R(:,i-delka_okna)*1/delka_okna;
    end 
    figure(1)
    stem(kum)
    xlim([0,4])
    ylim([0,10])
    pause(1)
end 

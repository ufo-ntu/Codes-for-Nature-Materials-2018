clear 
clc

fs = 60;           %sampling rate
ts = 1/fs;
np9 = 35;
tune = 1;
T = 0.002;
ratio_peak = 0.8893;
necho = [1 370 295 1035];                          


new = dlmread(' ');  %file name for A1
new1 = dlmread(' '); %file name for A3
time = new(:,1);
signal = new(:,2);

b0 = [0.4 1.43];
b = fminsearch(@simulation_fit_sample2,b0,[],signal,time,np9,tune);              %nonlinear fitting

period = 1/b(2);
t = 0:ts:20;
duty = 2.4/(2.4+3.1)*100;
sens = (square(t*2*pi/period,duty) + 1)/2;
se = find(abs(diff(sens)) == 1);
sens = sens(:,[se(1)+1:se(20)]);

test = (square(t*2*pi/period) + 1)/2;
wa = find(abs(diff(test)) == 1);
wave = cos(2*pi/period.*t);
wave = wave(:,[wa(1)+1:wa(19)+1]);
sim3 = conv(wave,sens);
sim1 = sim3*2;
valley = find(sim1(1:end-2)>sim1(2:end-1) & sim1(2:end-1)<sim1(3:end)) + 1;
peak = find(sim1(1:end-2)<sim1(2:end-1) & sim1(2:end-1)>sim1(3:end)) + 1;
sim1(:,[1:peak(9)-1]) = [];   
sim1 = sim1.*10^-9*b(1);
O = zeros(1,np9-1);
sim1 = [O sim1];
sim1([1,length(sim1)+1:length(signal)]) =  0;
sim1 = sim1';

bg = signal - sim1;
top = find(bg == max(bg))+tune;
bg = bg([top:end],:);                                 
sig11 = signal([top:end],:); 
t11 = time([top:end],:);

s1 = bg.*10^6;
c = polyfit(t11,s1,9);
fit1 = c(1).*t11.^9+c(2).*t11.^8+c(3).*t11.^7+c(4).*t11.^6+c(5).*t11.^5+c(6).*t11.^4+c(7).*t11.^3+c(8).*t11.^2+c(9).*t11+c(10);
fit1 = fit1./10^6;
bp1 = sig11 - fit1;
bg1 = bg - fit1;

plot(t11,bp1)

a = 4;
decay1 = exp(-a*(t11-t11(necho(2)))).*heaviside(t11 - t11(necho(2)))+exp(a*(t11-t11(necho(1)))).*-(heaviside(t11 - t11(necho(1)))-1);
decay1([necho(1):necho(2)],:) = 1;
echo1 = bp1.*decay1;
long = fs*1000;
f = (0:long-1)*fs/long;
echo1([length(echo1)+1:long],:) =  0;
bg1([length(bg1)+1:long],:) =  0;
f1 = abs(fft(echo1));
fbg1 = abs(fft(bg1));
n1 = find(abs(f-1.2) == min(abs(f-1.2)));
n2 = find(abs(f-1.6) == min(abs(f-1.6)));

error = (sum(f1([n1:n2],1)))*1e3;

[b(1) b(2) max(f1)*1e6] 

%%

cf = b(2);
b = b(1);
first = 0;


%starin simulation
period = 1/cf;
fs = 60;
ts = 1/fs;
long = fs*1000;
t = 0:ts:20;
duty = 2.4/(2.4+3.1)*100;
sens = (square(t*2*pi/period,duty) + 1)/2;
se = find(abs(diff(sens)) == 1);
sens = sens(:,[se(1)+1:se(20)]);

test = (square(t*2*pi/period) + 1)/2;

wa = find(abs(diff(test)) == 1);
wave = cos(2*pi/period.*t);
wave = wave(:,[wa(1)+1:wa(19)+1]);
sim3 = conv(wave,sens);
sim1 = sim3*2;
valley = find(sim1(1:end-2)>sim1(2:end-1) & sim1(2:end-1)<sim1(3:end)) + 1;
peak = find(sim1(1:end-2)<sim1(2:end-1) & sim1(2:end-1)>sim1(3:end)) + 1;
sim1(:,[1:valley(10)+first]) = [];   

f = (0:long-1)*fs/long;
sim1(:,[length(sim1)+1:long]) =  0;
sim3(:,[length(sim3)+1:long]) =  0;
ss = sim3*2;
ss(:,[1:peak(9)-1]) = [];
O = zeros(1,np9-1);
ss = [O ss];



fs = 60;           %sampling rate
ts = 1/fs;
scan = ceil(0.2/ts);

%load A1
time = new(:,1);
sig = new(:,2);
B1 = max(sig);
ss = ss(:,[1:length(sig)])'*10^-9*b;
cd = sig - ss;                                              
top = find(cd == max(cd))+tune;
cd = cd([top:end],:);
sig11 = sig([top:end],:);                                 
t11 = time([top:end],:);

test1 = cd;
s1 = test1.*10^6;
c = polyfit(t11,s1,9);
fit1 = c(1).*t11.^9+c(2).*t11.^8+c(3).*t11.^7+c(4).*t11.^6+c(5).*t11.^5+c(6).*t11.^4+c(7).*t11.^3+c(8).*t11.^2+c(9).*t11+c(10);
fit1 = fit1./10^6;
bp1 = sig11 - fit1;

sm1 = fit1;
f1 = (0:length(sig11)-1)*fs/length(sig11);
fsig1 = fft(sig11);
fsm1 = fft(sm1);


%load A3
t13 = new1(:,1);
sig13 = new1(:,2);
n = find(sig13 == max(sig13));
t13 = t13([n:end],:);
sig13 = sig13([n:end],:);


test3 = sig13;
sm3 = smooth(test3,scan);
for i = 1:1:100
    sm3 = smooth(sm3,scan);
end
f3 = (0:length(sig13)-1)*fs/length(sig13);
ftest3 = fft(test3); 
fsm3 = fft(sm3);
bp3 = test3 - sm3;


figure(1)
subplot(2,2,1);plot(t11,sig11,t11,fit1);
subplot(2,2,2);plot(bp1);
subplot(2,2,3);plot(f1,abs(fsig1),f1,abs(fsm1));axis([0 2 0 1e-3]);
subplot(2,2,4);plot(t11,test1,t11,fit1);

figure(2)
subplot(2,2,1);plot(t13,test3,t13,sm3);
subplot(2,2,2);plot(bp3);
subplot(2,2,3);plot(f3,abs(ftest3),f3,abs(fsm3));axis([0 2 0 1e-4]);
subplot(2,2,4);plot(t13,bp3);

bp1 = bp1./T;
bp3 = bp3./T;
sim1 = sim1./T;

%% lifetime analyze
a = 4;
decay1 = exp(-a*(t11-t11(necho(2)))).*heaviside(t11 - t11(necho(2)))+exp(a*(t11-t11(necho(1)))).*-(heaviside(t11 - t11(necho(1)))-1);
decay1([necho(1):necho(2)],:) = 1;
decay3 = exp(-a*(t13-t13(necho(4)))).*heaviside(t13 - t13(necho(4)))+exp(a*(t13-t13(necho(3)))).*-(heaviside(t13 - t13(necho(3)))-1);
decay3([necho(3):necho(4)],:) = 1;
echo1 = bp1.*decay1;
echo3 = bp3.*decay3;

sim1 = sim1(:,[1:length(sig11)])'*10^-9*b;
figure(4)
subplot(2,2,1);plot(t11,echo1);
subplot(2,2,2);plot(t11,echo1,t11,sim1);
subplot(2,2,3);plot(t13,echo3);

long = fs*1000;
f = (0:long-1)*fs/long;
echo1([length(echo1)+1:long],:) =  0;
echo3([length(echo3)+1:long],:) =  0;
sim1([length(sim1)+1:long],:) =  0;
f1 = fft(echo1);
f3 = fft(echo3);
fsim1 = fft(sim1);
subplot(2,2,4);plot(f,abs(f1),f,abs(fsim1));
axis([0 2 0 15e-3]);


nmax1 = find(abs(f1)==max(abs(f1)));
nmax1 = nmax1(1);
left1 = abs(abs(f1([1:nmax1],:)) - max(abs(f1))/2);
right1 = abs(abs(f1([nmax1:2500],:)) - max(abs(f1))/2);
w1 = [find(left1 == min(left1)) find(right1 == min(right1))+nmax1-1];
F1 = abs(f1(w1(1):w1(2),1));

nmax3 = find(abs(f3)==max(abs(f3)));
nmax3 = nmax3(1);
left3 = abs(abs(f3([1:nmax3],:)) - max(abs(f3))/2);
right3 = abs(abs(f3([nmax3:2500],:)) - max(abs(f3))/2);
w3 = [find(left3 == min(left3)) find(right3 == min(right3))+nmax3-1];
F3 = abs(f3(w3(1):w3(2),1));

T13 = 141.4;

f1_rv1_pv1_f3_rv3_pv3 = [f(find(F1 == max(F1))+w1(1)-1),sum(F1)/length(F1)*1e3,max(F1)*1e3,f(find(F3 == max(F3))+w3(1)-1),sum(F3)/length(F3)*1e3,max(F3)*1e3]
peak_ratio_tau = [(max(F3))/(max(F1))  T13./(2*log(max(F1)/max(F3)/ratio_peak))]


figure(5)
subplot(2,2,1);plot(f,abs(f1),'k',f,abs(f3),'r');axis([0 2 0 5e-3]);
subplot(2,2,3);plot(f,abs(f1),'k');
line([f(w1(1)) f(w1(1))],ylim);line([f(w1(2)) f(w1(2))],ylim);axis([0 2 0 5e-3]);
subplot(2,2,4);plot(f,abs(f3),'r');
line([f(w3(1)) f(w3(1))],ylim);line([f(w3(2)) f(w3(2))],ylim);axis([0 2 0 5e-3]);

%%
%data output normalized
n = find(f == 3);
f = f(:,[1:n]);
F1 = abs(f1([1:n],:));
norm = max(F1);
F1 = F1./norm;
F3 = abs(f3([1:n],:));
F3 = F3./norm;
echo = [f;F1';f;F3'];
fid = fopen('A1570_freq2.txt','w');                             
fprintf(fid,'%9f  %12f  %15f  %18f\r\n',echo);
fclose(fid);

A1 = [t11';bp1'.*10^3./norm];
A3 = [t13';bp3'.*10^3./norm];

fid = fopen('A1570_A1.txt','w');                              
fprintf(fid,'%9f  %12f\r\n',A1);
fclose(fid);

fid = fopen('A1570_A3.txt','w');                              
fprintf(fid,'%9f  %12f\r\n',A3);
fclose(fid);

%data output_original 

data1 = [t11';bp1'.*10^5];
data3 = [t13';bp3'.*10^5];

fid = fopen('A1570_data11.txt','w');                              
fprintf(fid,'%9f  %12f\r\n',data1);
fclose(fid);

fid = fopen('A1570_data13.txt','w');                             
fprintf(fid,'%9f  %12f\r\n',data3);
fclose(fid);

%%
time = new(:,1);
sig = new(:,2)./T;
trace = [time';sig'.*10^5];
fid = fopen('A1570_trace11.txt','w');                              
fprintf(fid,'%9f  %12f\r\n',trace);
fclose(fid);

time = new1(:,1);
sig = new1(:,2)./T;
trace = [time';sig'.*10^5];
fid = fopen('A1570_trace13.txt','w');                              
fprintf(fid,'%9f  %12f\r\n',trace);
fclose(fid);
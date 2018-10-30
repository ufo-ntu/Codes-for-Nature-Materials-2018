clear
clc

b = 1;
first = 1;
cf = 1.06;

period = 1/cf;
fs = 60;
ts = 1/fs;
long = fs*1000;
t = 0:ts:20;
duty = 3/(8/cf)*100;
sens = (square(t*2*pi/period,duty) + 1)/2;
se = find(abs(diff(sens)) == 1);
sens = sens(:,[se(1)+1:se(22)]);

test = (square(t*2*pi/period) + 1)/2;
wa = find(abs(diff(test)) == 1);
wave = cos(2*pi/period.*t);
wave = wave(:,[wa(1)+1:wa(21)]);

sim3 = conv(wave,sens);
sim1 = sim3*2;
valley = find(sim1(1:end-2)>sim1(2:end-1) & sim1(2:end-1)<sim1(3:end)) + 1;
peak = find(sim1(1:end-2)<sim1(2:end-1) & sim1(2:end-1)>sim1(3:end)) + 1;
sim1o = sim1;
sim1(:,[1:valley(11)+first]) = []; 
sim1o(:,[1:peak(10)-1]) = []; 

figure(2)
subplot(2,2,1);plot(sim1);
subplot(2,2,2);plot(sim3);
subplot(2,2,3);plot(sim1o);
%%

f = (0:long-1)*fs/long;
sim1(:,[length(sim1)+1:long]) =  0;
sim1o(:,[length(sim1o)+1:long]) =  0;
sim3(:,[length(sim3)+1:long]) =  0;
f1 = fft(sim1);
f1o = fft(sim1o);
f3 = fft(sim3);

f1 = f1';
f1o = f1o';
f3 = f3';

figure(5)
subplot(2,2,1);plot(f,abs(f1),'k',f,abs(f3),'r');axis([0 2 0 5e4]);
subplot(2,2,2);plot(f,abs(f1o),'k',f,abs(f3),'r');axis([0 2 0 5e4]);

n = find(f == 3);
f = f(:,[1:n]);
F1 = abs(f1([1:n],:)).*10^-4;
F1o = abs(f1o([1:n],:)).*10^-4;
F3 = abs(f3([1:n],:)).*10^-4;

ratio_peak = max(F1)/max(F3)
ratio_peak_o = max(F1o)/max(F3)

echo = [f;F1';f;F3'];
echoo = [f;F1o';f;F3'];
fid = fopen('sim_freq.txt','w');                              %輸出檔案的名稱##########################################
fprintf(fid,'%9f  %12f  %15f  %18f\r\n',echo);
fclose(fid);
fid = fopen('simo_freq.txt','w');                              %輸出檔案的名稱##########################################
fprintf(fid,'%9f  %12f  %15f  %18f\r\n',echoo);
fclose(fid);
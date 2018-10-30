function F = strain_fit(b,signal,time,np9,tune)
fs = 60;           %sampling rate
ts = 1/fs;

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
sim1(:,[1:peak(9)-1]) = [];   %8
sim1 = sim1.*10^-9*b(1);
O = zeros(1,np9-1);
sim1 = [O sim1];
sim1([1,length(sim1)+1:length(signal)]) =  0;
sim1 = sim1';

bg = signal - sim1;
top = find(bg == max(bg))+tune;
bg = bg([top:end],:);                                 
sig1 = signal([top:end],:); 
t1 = time([top:end],:);

s1 = bg.*10^6;
c = polyfit(t1,s1,9);
fit1 = c(1).*t1.^9+c(2).*t1.^8+c(3).*t1.^7+c(4).*t1.^6+c(5).*t1.^5+c(6).*t1.^4+c(7).*t1.^3+c(8).*t1.^2+c(9).*t1+c(10);
fit1 = fit1./10^6;
bp1 = sig1 - fit1;
bg1 = bg - fit1;


necho = [1 320];                          
a = 4;
decay1 = exp(-a*(t1-t1(necho(2)))).*heaviside(t1 - t1(necho(2)))+exp(a*(t1-t1(necho(1)))).*-(heaviside(t1 - t1(necho(1)))-1);
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

F = (sum(fbg1([n1:n2],1))-sum(f1([n1:n2],1)))*1e3;
clear 
[x, fsample]= audioread('musicex1.wav');

x = x(fsample*1:fsample*10)
N = length(x);
N1 = 1000;
N2 = 6;

delta_f = fsample/N;
f_axis = [0:delta_f:fsample-delta_f];

%Lavpas
fc_LP = 4000;

%Båndpas1
fc_BP1 =3900;
fc_BP11= 8800;

%Båndpas2 
fc_BP2 = 8700;
fc_BP22 = 13200;

%Båndpas3 
fc_BP3 = 13100;
fc_BP33 = 17600;

%Højpas
fc_HP = 17500;


k1 = 1;
k2 = 1;
k3 = 50;
k4 = 0.1;
k5 = 1;

%%
%Filtre

%Lavpas
%LP = fir1(N1, fc_LP/(0.5*fsample));
[b,a] = butter(N2,fc_LP/(0.5*fsample));

%Højpas
HP = fir1(N1, fc_HP/(0.5*fsample), 'high');

%Båndpas
BP1 = fir1(N1, [fc_BP1 fc_BP11]/(0.5*fsample), 'bandpass');
BP2 = fir1(N1, [fc_BP2 fc_BP22]/(0.5*fsample), 'bandpass');
BP3 = fir1(N1, [fc_BP3 fc_BP33]/(0.5*fsample), 'bandpass');

%%
%Filtrering
y_LP = filter(b,a,x);
y_HP = filter(HP,1,x);
y_BP1 = filter(BP1,1,x);
y_BP2 = filter(BP2,1,x);
y_BP3 = filter(BP3,1,x);


%y_LP=conv(x,LP);
%y_HP=conv(x,HP);
%y_BP1=conv(x,BP1);
%y_BP2=conv(x,BP2);
%y_BP3=conv(x,BP3);

%%
%Samlede output for signalet
y_EQ = k1*y_LP + k2*y_HP + k3*y_BP1 + k4*y_BP2 + k5*y_BP3;

%Her kan man se, hvordan de forskellige filtre overlapper hinanden. 
%Det skal gerne være en lige streg, som man kan se nu.

figure(1)
freqz(b,a)
hold on
freqz(BP1)
freqz(BP2)
freqz(BP3)
freqz(HP)

figure(2)
freqz(y_EQ)
title('Det filtreret signal')

%%
%FFT lavpas
X = fft(x, N);
[h,t] = impz(b,a,50);

figure(17)
plot(h)
xlabel('Tid')
ylabel('Amplitude')

LPf = fft(h.', N);
Y_LPf = X.*LPf;

%%
figure(3)
semilogx(f_axis(1:0.5*end),20*log10(abs(2/N)*Y_LPf(1:0.5*end)))
title('LP FFT')
xlabel('Hz')
ylabel('dB')

figure(4)
plot(f_axis(1:0.5*end), unwrap(angle(Y_LPf(1:0.5*end)))) %Unwrap bruges til at lægge faserne sammen
xlabel('Hz')
ylabel('Fase')

%FFT båndpas1
BP1f = fft(BP1, N);
Y_BP1f = X.*BP1f;

figure(5)
semilogx(f_axis(1:0.5*end),20*log10(abs(2/N)*Y_BP1f(1:0.5*end)))
title('Bp1 FFT')
xlabel('Hz')
ylabel('dB')

figure(6)
plot(f_axis(1:0.5*end), unwrap(angle(Y_BP1f(1:0.5*end)))) %Unwrap bruges til at lægge faserne sammen
xlabel('Hz')
ylabel('Fase')

%FFT båndpas2
BP2f = fft(BP2, N);
Y_BP2f = X.*BP2f;

figure(7)
semilogx(f_axis(1:0.5*end),20*log10(abs(2/N)*Y_BP2f(1:0.5*end)))
title('Bp2 FFT')
xlabel('Hz')
ylabel('dB')

figure(8)
plot(f_axis(1:0.5*end), unwrap(angle(Y_BP2f(1:0.5*end)))) %Unwrap bruges til at lægge faserne sammen
xlabel('Hz')
ylabel('Fase')

%FFT båndpas3
BP3f = fft(BP3, N);
Y_BP3f = X.*BP3f;

figure(9)
semilogx(f_axis(1:0.5*end),20*log10(abs(2/N)*Y_BP3f(1:0.5*end)))
title('Bp3 FFT')
xlabel('Hz')
ylabel('dB')

figure(10)
plot(f_axis(1:0.5*end), unwrap(angle(Y_BP3f(1:0.5*end)))) %Unwrap bruges til at lægge faserne sammen
xlabel('Hz')
ylabel('Fase')

%FFT højpas
HPf = fft(HP, N);
Y_HPf = X.*HPf;

figure(11)
semilogx(f_axis(1:0.5*end),20*log10(abs(2/N)*Y_HPf(1:0.5*end)))
title('Hp FFT')
xlabel('Hz')
ylabel('dB')

figure(12)
plot(f_axis(1:0.5*end), unwrap(angle(Y_HPf(1:0.5*end)))) %Unwrap bruges til at lægge faserne sammen
xlabel('Hz')
ylabel('Fase')

%FFT samlede output 
Y_EQf = k1*Y_LPf + k2*Y_HPf + k3*Y_BP1f + k4*Y_BP2f + k5*Y_BP3f;

figure(13)
semilogx(f_axis(1:0.5*end),20*log10(abs(2/N)*Y_EQf(1:0.5*end)))
title('FFT for siganlet, igennem alle filtre')
xlabel('Hz')
ylabel('dB')

figure(14)
plot(f_axis(1:0.5*end), unwrap(angle(Y_EQf(1:0.5*end)))) %Unwrap bruges til at lægge faserne sammen
xlabel('Hz')
ylabel('Fase')

figure(16)
freqz(x)
title('Ufiltreret')

%%
%Invers af fft

%Lavpas
EQi = ifft(Y_EQf);

figure(15)
freqz(EQi)
title('Inverse FFT')
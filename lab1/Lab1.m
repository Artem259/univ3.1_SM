%#ok<*NOPTS>
clear % видалення попередніх змінних

T = 5; % період сигналу
Fs = 1/T; % частота сигналу             
dT = 0.01;

Y = dlmread('f8.txt', ' '); % читання файлу у масив вимірів сигналу
Y = Y(1:end-1);
L = length(Y);

t = (0:L-1)*dT; % формування масиву моментів часу

tf = (0:L-1); % формування масиву частот
tf = tf(1:L/2);

figure('Name','Data initial plot');
plot(t, Y), grid

Yf = fft(Y)/L; % перетворення Фур'є
Yf = abs(Yf);
Yf = Yf(1:L/2);

figure('Name','Fourier transform plot');
plot(tf,Yf), grid

% локальні максимуми модуля перетворення Фур'є:
i = islocalmax(Yf);
i = tf(i);
% ---

f = i*Fs  

f_sin = sin(2*pi*f*t);

M = [sum(t.^6),        sum(t.^5),        sum(t.^4),     sum(f_sin.*t.^3),  sum(t.^3);
     sum(t.^5),        sum(t.^4),        sum(t.^3),     sum(f_sin.*t.^2),  sum(t.^2);
     sum(t.^4),        sum(t.^3),        sum(t.^2),     sum(f_sin.*t),     sum(t);
     sum(f_sin.*t.^3), sum(f_sin.*t.^2), sum(f_sin.*t), sum(f_sin.*f_sin), sum(L*f_sin);
     sum(t.^3),        sum(t.^2),        sum(t),        sum(L*f_sin),      L];

c = [sum(Y.*t.^3), sum(Y.*t.^2), sum(Y.*t), sum(Y.*f_sin),  sum(Y)];

A = (M\c')'

resY = A(1).*t.^3 + A(2).*t.^2 + A(3).*t + A(4).*f_sin + A(5);

F = sum((resY-Y).^2)/2 % функціонал похибки

figure('Name','Data approximate plot');
plot(t, resY), grid

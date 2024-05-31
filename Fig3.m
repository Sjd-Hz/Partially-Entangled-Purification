close all
clear
clc
a=linspace(0,sqrt(2)/2,100);
y1=4*a.^2.*(1-a.^2);
y2=2*a.^2;
plot(y1,y2)
xlabel("Linear Entropy S(\psi)")
ylabel("Max(P_{sum}^{net})")
ylim([0 1])
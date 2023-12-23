clear all
close all
clc
%% Hopf bifurcation curve 
p = 0.2; %default parameters for WC system 
a = 5;
tau1 = 1;
tau2 = 5;
alpha = 1;

% curve was analytically derived from Nicola et. al., 2018, Chaos.  See
% citation for furhter details.
theta = 0:0.001:10;
thetap = 1./(1+exp(-a*theta*p));

D = (thetap.^2)*a*p*(1-p)/(tau1*tau2); 
K = p*a*thetap.*(1-thetap).*theta./(thetap);
pinv = (1/a)*(log(p/(1-p))); 
F = (1/alpha)/tau1 -  ((log(p/(1-p)))/a)*a*p*(1-p)/(p*tau1);

b1 = (K.*F+D+1-K/alpha);
a1 = (1-K/alpha); 
c1 = K.*F;
mup = (-b1 + sqrt(b1.^2-4*a1.*c1))./(2*a1);
mum = (-b1 -sqrt(b1.^2-4*a1.*c1))./(2*a1);
WEp = (1-tau1*mup)/(a*p*(1-p)*alpha);
WEm = (1-tau1*mum)/(a*p*(1-p)*alpha);
subplot(4,2,[1,3,5,7])
plot(theta,WEp,'LineWidth',2), hold on
ylim([1.1,2.3])
xlabel('$W^{IE}$','Interpreter','LaTEx','FontSize',14)
xlabel('$W^{E}$','Interpreter','LaTEx','FontSize',14)


%% Pre-Hopf simulation and plotting 
figure(1)
options = odeset('AbsTol',1e-14,'RelTol',1e-14);
theta =1;
WE = 1.9;
T = 10^3;
[tsn2,ysn2]= ode45(@(t,y) WC(WE,tau1,tau2,p,a,theta,1,y),0:0.1:T,rand(3,1),options);
[t,y]= ode45(@(t,y) WC(WE,tau1,tau2,p,a,theta,1,y),0:0.1:T,ysn2(end,:)',options);
subplot(4,2,2)
plot(t,y(:,1),'LineWidth',2)
ylim([0,1])
subplot(4,2,[1,3,5,7])
plot(theta,WE,'r.','MarkerSize',14)
xlabel('$W^{IE}$','Interpreter','LaTEx','FontSize',14)
ylabel('$W^{E}$','Interpreter','LaTEx','FontSize',14)
figure(2)
subplot(1,4,1)
plot(y(:,1),y(:,2)) 
xlim([0,1])
ylim([0,1])
xlabel('$E(t)$','Interpreter','LaTEx','FontSize',14)
ylabel('$I(t)$','Interpreter','LaTEx','FontSize',14)
title('$W^E = 1.9$','Interpreter','LaTex','FontSize',14)

%% Post-Hopf 
figure(1)
options = odeset('AbsTol',1e-14,'RelTol',1e-14);
theta =1;
WE = 2.05;
T = 10^3;
[tsn2,ysn2]= ode45(@(t,y) WC(WE,tau1,tau2,p,a,theta,1,y),0:0.1:T,rand(3,1),options);
[t,y]= ode45(@(t,y) WC(WE,tau1,tau2,p,a,theta,1,y),0:0.1:T,ysn2(end,:)',options);
subplot(4,2,4)
plot(t,y(:,1),'LineWidth',2)
xlabel('$t$','Interpreter','LaTEx','FontSize',14)
ylabel('$E(t)$','Interpreter','LaTEx','FontSize',14)

ylim([0,1])
subplot(4,2,[1,3,5,7])
plot(theta,WE,'r.','MarkerSize',14)
figure(2)
subplot(1,4,2)
plot(y(:,1),y(:,2)) 
xlim([0,1])
ylim([0,1])
xlabel('$E(t)$','Interpreter','LaTEx','FontSize',14)
ylabel('$I(t)$','Interpreter','LaTEx','FontSize',14)

title('$W^E = 2.05$','Interpreter','LaTex','FontSize',14)

%% Chaotic System 
figure(1)
options = odeset('AbsTol',1e-14,'RelTol',1e-14);
theta =1;
WE = 2.115;
T = 10^3;
[tsn2,ysn2]= ode45(@(t,y) WC(WE,tau1,tau2,p,a,theta,1,y),0:0.1:T,rand(3,1),options);
[t,y]= ode45(@(t,y) WC(WE,tau1,tau2,p,a,theta,1,y),0:0.1:T,ysn2(end,:)',options);
subplot(4,2,6)
plot(t,y(:,1),'LineWidth',2)
ylim([0,1])
subplot(4,2,[1,3,5,7])
plot(theta,WE,'r.','MarkerSize',14)
figure(2)
subplot(1,4,3)
plot(y(:,1),y(:,2)) 
xlim([0,1])
ylim([0,1])
xlabel('$E(t)$','Interpreter','LaTEx','FontSize',14)
ylabel('$I(t)$','Interpreter','LaTEx','FontSize',14)

title('$W^E = 2.115$','Interpreter','LaTex','FontSize',14)

%% Relaxation Oscillator
figure(1)
options = odeset('AbsTol',1e-14,'RelTol',1e-14);
theta =1;
WE = 2.25;
T = 10^3;
[tsn2,ysn2]= ode45(@(t,y) WC(WE,tau1,tau2,p,a,theta,1,y),0:0.1:T,rand(3,1),options);
[t,y]= ode45(@(t,y) WC(WE,tau1,tau2,p,a,theta,1,y),0:0.1:T,ysn2(end,:)',options);
subplot(4,2,8)
plot(t,y(:,1),'LineWidth',2)
ylim([0,1])
subplot(4,2,[1,3,5,7])
plot(theta,WE,'r.','MarkerSize',14)
figure(2)
subplot(1,4,4)
plot(y(:,1),y(:,2)) 
xlim([0,1])
ylim([0,1])
xlabel('$E(t)$','Interpreter','LaTEx','FontSize',14)
ylabel('$I(t)$','Interpreter','LaTEx','FontSize',14)

title('$W^E = 2.25$','Interpreter','LaTex','FontSize',14)


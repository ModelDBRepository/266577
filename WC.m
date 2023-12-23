function dy = WC(WE,tau1,tau2,p,a,theta,N,y)
E=y(1:N); 
I=y(N+1:2*N); 
W=y(2*N+1:3*N); 
A1 = WE*E - W.*I;
A2 = theta*E; 
dE = (-E + 1./(1+exp(-a*(A1))))/tau1;
dI = (-I + 1./(1+exp(-a*A2)));
dW = I.*(E-p)/tau2;
dy = [dE;dI;dW];
end


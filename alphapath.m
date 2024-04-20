%EM Euler-Maruyama method on linear USDE
%
% USDE is dX = lambda*X dt + mu*X dW + , X(0) = Xzero,
% Discretized Brownian path over [0,1] has dt = 2^(-8).
% Euler-Maruyama uses timestep R*dt.
r=rng('default',100);
lambda = 0.06; sigma1 = 0.32; sigma2 = 0.29; afa = 0.95; Xzero = 40; % problem parameters
T = 1; N = 100; dt = 1/N; alpha = ((2*sqrt(3)/pi)*log(afa/(1-afa)));
dW = sqrt(dt)*rng(1,N); % Brownian increments
W = cumsum(dW); % discretized Brownian path


R = 1; Dt = R*dt; L = N/R; % L EM steps of size Dt = R*dt
Xem = zeros(1,L); % preallocate for efficiency
Xtemp = Xzero;
for j = 1:L
Winc = sum(dW(R*(j-1)+1:R*j));
Xtemp = Xtemp + Dt*lambda*Xtemp + sigma1*Xtemp*Winc + Dt*abs(sigma2*Xtemp)*alpha;
Xem(j) = Xtemp;
end
plot([0:Dt:T],[Xzero,Xem],'r--*'), hold off
xlabel('t','FontSize',12)
ylabel('X','FontSize',16,'Rotation',0,'HorizontalAlignment','right')
expectedUSDE = sum(Xem)/N
%EM Euler-Maruyama method on linear USDE
%
% USDE is dX = lambda*X dt + mu*X dW, X(0) = Xzero,
% where lambda = 2, mu = 1 and Xzero = 1.
%
% Discretized Brownian path over [0,1] has dt = 2^(-8).
% Euler-Maruyama uses timestep R*dt.
randn('state',100)
lambda = 0.06; mu = 0.32; v = 0.29; afa = 0.50; Xzero = 40; % problem parameters
T = 1; N = 100; dt = 1/N; alpha = ((2*sqrt(3)/pi)*log(afa/(1-afa)));
dW = sqrt(dt)*randn(1,N); % Brownian increments
W = cumsum(dW); % discretized Brownian path
Xtrue = Xzero*exp((lambda-0.5*mu^2)*([dt:dt:T])+(mu*W) + v*alpha);%Exact solution
plot([0:dt:T],[Xzero,Xtrue],'m-'), hold on
R = 1; Dt = R*dt; L = N/R; % L EM steps of size Dt = R*dt
Xem = zeros(1,L); % preallocate for efficiency
Xtemp = Xzero;
for j = 1:L
Winc = sum(dW(R*(j-1)+1:R*j));
Xtemp = Xtemp + Dt*lambda*Xtemp + mu*Xtemp*Winc + Dt*abs(v*Xtemp)*alpha;
Xem(j) = Xtemp;
end
plot([0:Dt:T],[Xzero,Xem],'r--*'), hold off
xlabel('t','FontSize',12)
ylabel('X','FontSize',16,'Rotation',0,'HorizontalAlignment','right')
emerr = abs(Xem(end)-Xtrue(end))
expectedXact = sum(Xtrue)/N
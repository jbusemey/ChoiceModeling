function [pr, M, pdf, cdf] = DFT2ALTCRT(parms,T)
% function [pr, M, pdf, cdf] = DFT2ALTCRT(parms,T)
% Choice and Response Time Dist for 2 Alt DFT 
% Busemeyer & Townsend (1993, Psychological Review)

% comment out line 118 to compute cdf if not needed to speed up calculation 

% Diffusion model parms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mu = mean drift rate
% std = sqrt( diffusion rate)
% d = mu/std   (only need to input this ratio)
%    positive value drives process down matrix
% theta = threshold (in std units)
%     evidence ranges from -theta (top) to +theta (bottom) of matrix
% z = bias (in std units)
%    positive value bias process down matrix
% s = growth-decay (s = 0 is Weiner, otherwise OU)

% input  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parms = [d theta z s];
% example parameters
%           d  theta   bias    s
% parms = [.5    1     -.2     0];

% T = N x 1 vector of decision times (designed for seconds)


% output
% up is alternative chosen at -theta bound top    of matrix
% dn is alternative chosen at +theta bound bottom of matrix

% pr = choice prob for each alt (up alt, dn alt)
% M = mean choice time conditioned on each alt (up alt, dn alt)
% pdf = prob density stop at time t and choose for each alt (up,dn)
%      vector with N rows representing the N input times
% cdf = cum distribution at time t for each alt (up,dn)
%      vector with N rows representing the N input times
%      cdf(up at t=inf) = pr(up), cdf(dn at t=inf) = pr(dn)
%      cdf(up)+cdf(dn) = 1

d = parms(1);        %  .5; pos -> drift down matrix towards increasing values
theta = parms(2);    %  1;
z = parms(3);        %   -.20;  neg -> bias up toward neg values above zero
s = parms(4);        %    0;  % Weiner process

% build start   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = .0004;       % time step
del = sqrt(h);   % step size
Dist = theta/del;
Mid = ceil(Dist);  % no steps to threshold
ns = 2*Mid +1;
X = (-Mid:Mid)'.*del;

S0 = zeros(ns,1);
% smooth initial bias to avoid discreteness
bias = z/del;
bias = Mid+1+bias;
fb = floor(bias);
cb = ceil(bias);
if fb ~= cb
    S0(fb) = bias-fb;
    S0(cb) = cb-bias;
else
    S0(cb) = 1;
end
S0 = S0(2:(ns-1),1);

% Continuous time discrete state Markov model %%%%%%%%%%%%%%%%%%%%%%%%%%%

dX = d - s*X;
a = 0;
q = (1./(2+a)).*(1 + dX.*del);  % dn prob
p = (1./(2+a)).*(1 - dX.*del);  % up prob
r = 1 - (p+q);
if (sum(q<0)||sum(q>1)||sum(p<0)||sum(p>1)) 
   'invalid s parameter'
end
k = 1/del;

gam = -k^2;
alfa = k^2.*p./(p+q);    % rate up = (prob up)/h
beta = k^2.*q./(p+q);    % rate dn = (prob dn)/h

up = alfa(2:ns,1)  ;
dn = beta(1:(ns-1),1);
cn = gam.*ones(ns,1);

A = diag(up);
B = diag(dn);
Q = diag(cn);
Q(1:(ns-1),2:ns) = Q(1:(ns-1),2:ns) + A;     
Q(2:ns,1:(ns-1)) = Q(2:ns,1:(ns-1)) + B;    % in col out row

% smooth threshold to avoid discreteness

fb = floor(Dist);
cb = ceil(Dist);
if fb ~= cb
   Q([1 2],3)          = [(cb-Dist) (Dist-fb)]'.*Q(2,3);
   Q([ns-1 ns],ns-2) =   [(Dist-fb) (cb-Dist)]'.*Q(ns-1,ns-2); 
end

% in col out row , col sums to zero
Q([1 2],1)= [ 0 0 ];
Q([ns-1 ns],ns) = [ 0 0];

Qs = Q(2:(ns-1),2:(ns-1));    %    in col out row  for DT, col sum = 0
Qt = Qs';                     %    in row out col  for pr, row sum = 0

% using to from row
Dq = diag(-diag(Qt).^-1);
K = Dq*Qt + eye(ns-2);          %  in row out col for pr, row sum = 1
R = zeros(ns-2,2);
R(1:2,1) = [Q(1,2) Q(1,3)]'/(-gam);        % choose up at -theta
R((ns-3):(ns-2),2) = [Q(ns,ns-2) Q(ns,ns-1)]'/(-gam);  % choose dn at +theta
PR = ((eye(ns-2)-K))\R;
pr = S0'*PR;
M = -PR'*(Qs\S0); M = M'./pr;

N = size(T,1);
pdf = zeros(N,2);
cdf = zeros(N,2);

for t=1:N
    % using to row from col
    pdf(t,:) = -PR'*Qs*expm(T(t).*Qs)*S0; 
    % comment out next line if cdf not needed
    cdf(t,:) = pr' -  PR'*expm(T(t).*Qs)*S0; 
    
end % t

end


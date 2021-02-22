function [sse, P3, TV] = fitMDFTs(x,D,FD,MM,C3,Ns)
% function [sse, P3, T] = fitMDFTs(x,D,FD,MM,C3,Ns)
% ccluates predicted and compare to observed

nc = size(D,1);
na = size(D,2);
P3 = zeros(nc,na);
TV = zeros(nc,1);

wgt = exp(x(1));
phi1 = exp(x(2)); 
phi2 = exp(x(3));
sig2 = exp(x(4));
theta1 = exp(x(5));
w = exp(x(6));
w = [w ; (1-w) ];  % weight vector for m attributes, m=2 in this case


for i = 1:nc
   
    M3 =MM(:,:,i);
    [G3, EG] = distfunct(M3,wgt,phi1,phi2);  % returns gamma
  %  EG
    [p3, T] = simMDF_mex(G3,C3,M3,w,theta1,sig2,Ns);   
   
    P3(i,:) =  p3;
    TV(i) = T;
  %   [P3 TV]
end

% disp( [wgt phi1 phi2 sig2 theta1 w' ])

 % choose to minimize sse or Chi
 
 % SSE
   dev = (D-P3);
   sse = sum(sum(dev.*dev));
   
%  Chi
%    
%    LL1 = FD.*log(P3); LL1 = sum(sum(LL1));
%    LL2 = FD.*log(D);  LL2 = sum(sum(LL2));
%    
%    Chi = -2*(LL2-LL1)
%    
   
  
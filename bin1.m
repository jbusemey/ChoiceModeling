function chi=bin1(parm,data)
% function chi=bin1(parm,data)
% computes chi square  for binomial model


a=parm(1); b = parm(2); c=parm(3); d = 1-(a+b+c);

y = data(:,1);
nt = size(y,1);

chi=0;

for t = 2:nt
yt = y(t,:) == [1 2 3 4]';
p = [a b c d]'; 

pp = .0001 + .9998*p;

chi = chi + log(pp)'*yt;
end % end of training

chi = -2*chi;

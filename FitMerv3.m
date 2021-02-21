function [ Chi ] = FitMerv3( Parm, PS, PM, nt )
% function [ Chi ] = FitMerv3( Parm, PS, PM, nt )
% three parm -- use only PD trials, bonus and not bonus trials

% W = 1./(1+exp(-Parm(1:5)));    % include for unconstrained

W = Parm(1);
a = Parm(2);
b = Parm(3); 
% b=0;
eps = 1./(10.^5);

Chi = 0;
B = PS(:,7);  % bonus indicator, 1 = not bonus, 2 = bonus
Op = PS(:,6); % Opp choice 

for row = 1:nt
    w1 = W; w2 = 1-w1; 
    w1 = (B(row)==1).*w1 + (B(row)==2).*(Op(row)==1);
    w2 = (B(row)==1).*w2 + (B(row)==2).*(Op(row)==2);
    Pay = PM(PS(row,3),:,PS(row,10));
    Pay = ( w1.*Pay(1) + w2.*Pay(3) ) - (w1.*Pay(5) + w2.*Pay(7));
    Pr = 1./(1+exp(-(a.*Pay+(B(row)==2).*b)));
    Pr = eps + (1-2.*eps).*Pr;
    ch = PS(row,4);   % obs choice
    LL = ch.*log(Pr) + (1-ch).*log(1-Pr);
    Chi = Chi + -2*LL;
end

end


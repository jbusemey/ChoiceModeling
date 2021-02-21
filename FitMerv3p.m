function [PP, Chi ] = FitMerv3p( Parm, PS, PM, nt )
% function [ Chi ] = FitMerv3p( Parm, PS, PM, nt )
% three parm -- use only PD trials, bonus and not bonus trials

% W = 1./(1+exp(-Parm(1:5)));    % include for unconstrained

W = Parm(1);
a = Parm(2);
b = Parm(3); 
% b=0;
eps = 1./(10.^5);

Chi = 0;
B = PS(:,7); 
Op = PS(:,6); 

% Pr1 = 0; Pr2 = 0; Ch1 = 0; Ch2 = 0;

PP = zeros(6,3);

for row = 1:nt
    w1 = W; w2 = 1-w1; 
    w1 = (B(row)==1).*w1 + (B(row)==2).*(Op(row)==1);
    w2 = (B(row)==1).*w2 + (B(row)==2).*(Op(row)==2);
    Pay = PM(PS(row,3),:,PS(row,10));
    Pay = ( w1.*Pay(1) + w2.*Pay(3) ) - (w1.*Pay(5) + w2.*Pay(7));
    Pr = 1./(1+exp(-(a.*Pay+(B(row)==2).*b)));
    Pr = eps + (1-2.*eps).*Pr;
    ch = PS(row,4);
    LL = ch.*log(Pr) + (1-ch).*log(1-Pr);
    Chi = Chi + -2*LL;
    
%     Pr1 = Pr1 + (B(row)==1).*Pr;    % no bonus
%     Pr2 = Pr2 + (B(row)==2).*Pr;    % bonus
%     Ch1 = Ch1 + (B(row)==1).*ch;
%     Ch2 = Ch2 + (B(row)==2).*ch;
     PP(PS(row,3),2) = PP(PS(row,3),2) + (B(row)==2).*(Op(row)==1).*Pr;
     PP(PS(row,3),3) = PP(PS(row,3),3) + (B(row)==2).*(Op(row)==2).*Pr;
     PP(PS(row,3),1) = PP(PS(row,3),1) + (B(row)==1).*Pr;
%       [PS(row,3) B(row) Op(row) ] 
%         PP(PS(row,3),:)

end

% Pr1 = Pr1/6;
% Pr2 = Pr2/12;
% Ch1 = Ch1/6;
% Ch2 = Ch2/12;
% Pr = [Pr1 Ch1 Pr2 Ch2];

end


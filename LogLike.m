% Compute log likelihoods
% build nsamples x nsub matrix  of log likelihoods

clear
clc
eps = 10^(-10);

load MervPD
% PD1 PM1

%  PD is N x 10 data matrix
%    1    2    3       4      5        6        7      8     9   10
%  subj trial matrix choice opponent oppchoice bonus payoff RT   Gnumb
% ns = 88; % no subj
% 18 trials per sub
% 12 bonus per person  bonus = 2
% 6  no bonus per person  bonus = 1
% 18*88 = 1584

% PM1  payoff matrices
% 6 rows x 8 col's matrix
% cols
% 1 2  3 4  5 6  7 8
% M Y  M Y  M Y  M Y
% D D  D C  C D  C C
% 1 1  1 2  2 1  2 2
% 1 = defect, 2 = coop

PD = reshape(PD1,18,88,10);
nsub = size(PD,2);
nt = size(PD,1);
PM = PM1;

    % change b to be nonzero
    load samplesI3_6  % 
    
%     rep = 4;
%     nr = 4;
%     ns = size(ParmsI,1)./nr;
%     nc = size(ParmsI,2);
%     
%     ParmsI = reshape(ParmsI,nr,ns,nc);
%     ParmsI = squeeze(ParmsI(rep,:,:));
%     
    
    
    % ParmsH
    % Wgt a b stda stdb : hyper means
    % ParmsI
    % wI aI bI   : for individuals
    
    %  1                 88
    % wI.1	wI.2 ...	wI.88    weight parm for subj i
    %  89                176
    % aI.1	...	        aI.88	 choice parm for subj i
    % 177                264
    % bI.1	bI.2...  	bI.88    bias parm for subj i
    
    ParmW = ParmsI(:,1:88);
    ParmA = ParmsI(:,89:176);
    ParmB = ParmsI(:,177:264);
    ns = size(ParmW,1);


%
Log_Like = zeros(ns,nsub);

for sub=1:nsub
    PS = squeeze(PD(:,sub,:));
    B = PS(:,7);
    Op = PS(:,6);
    
    for s = 1:ns
        w = ParmW(s,sub);
        a = ParmA(s,sub);
        b = ParmB(s,sub);
        
        
        LL=0;
        for row = 1:nt
            w1 = w; w2 = 1-w1;
            w1 = (B(row)==1).*w1 + (B(row)==2).*(Op(row)==1);
            w2 = (B(row)==1).*w2 + (B(row)==2).*(Op(row)==2);
            Pay = PM(PS(row,3),:,PS(row,10));
            Pay = ( w1.*Pay(1) + w2.*Pay(3) ) - (w1.*Pay(5) + w2.*Pay(7));
          % Pay = a.*( Pay + (B(row)==2).*b );   % I5
            Pay = a.*Pay + (B(row)==2).*b ;      % I6
            Pr = 1./(1+exp(-Pay));
            Pr = eps + (1-2.*eps).*Pr;
            ch = PS(row,4);
            LL = ch.*log(Pr) + (1-ch).*log(1-Pr);
        end
        Log_Like(s,sub) = LL;
    end  % end samples
end % end sub

 mstan.waic(Log_Like)
[loo3,loos3,pk3] = psisloo(Log_Like);

disp('Loo3')
disp(loo3)

% save Loo3 loos3

% sample I3- 5
% model 3 parm using not bestway to model b (a*b) model, sample 
% ans = 
% 
%          waic: 30.7125
%           lpd: -5.5175
%        p_waic: 9.8387
%     elpd_waic: -15.3562
%         p_loo: 10.4320
%      elpd_loo: -15.9496
% 
% Loo
%   -16.2556


% I3-6
% using added b  
% ans = 
% 
%          waic: 26.6827
%           lpd: -4.2011
%        p_waic: 9.1402
%     elpd_waic: -13.3413
%         p_loo: 9.4303
%      elpd_loo: -13.6315
% 
% Loo3
%   -13.6755


% ans = 
% 
%          waic: 26.7430
%           lpd: -4.2574
%        p_waic: 9.1141
%     elpd_waic: -13.3715
%         p_loo: 9.4092
%      elpd_loo: -13.6666
% 
% Loo3
%   -13.7317

% diff = loos3-loos2;
% sum(diff)=  3.6502
% std(diff)*sqrt(88)= 3.6292

 
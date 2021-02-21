% Stan Merv PD program
% 3 parameter model
% change how b is added 
% this is used for I6

clear
clc
load MervPD  % contains PD1 PM1

% PD1 data
% N x 10 data matrix
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

PM = PM1;
Nsub = 88;
N = size(PD1,1);
Nt = N/Nsub;

G = PD1(:,3);     % 1 to 6 matrices
PD = zeros(N,3);
PD(:,1) = 2 - PD1(:,4);   % player's: 1 = defect, 0 = coop
PD(:,2) = 2 - PD1(:,6);   % opp's:    1 = defect, 0 = coop
PD(:,3) = PD1(:,7)-1;     % trial type: 1 = bonus,  0 = no bonus

 LData = struct('N',N','Nsub',Nsub,'Nt',Nt,'PD',PD,'G',G,'PM',PM);
%  LData.N
 
MervModel = {
   'data {                                              '
   '    int<lower=0> Nsub; //No subj                    '
   '    int<lower=0> Nt; // number trials per sub       '
   '    int<lower=0> N ;   // Nsub*Nt                   '
   '    int<lower=0,upper=1> PD[N,3]; // Data           '
   '    int<lower=1,upper=6> G[N]; // Data              '
   '    int<lower=0,upper=30> PM[6,8]; // payoffs       '
   '}                                                   '
   
   'parameters {                                        '
   '    real<lower=1,upper=100>  mw1;                   '
   '    real<lower=1,upper=100>  mw2;                   '
   '    real ma ;                                       '
   '    real mb ;                                       '
   '    real<lower=0.01,upper=10> stda;                 '
   '    real<lower=0.01,upper=10> stdb;                 '
   '    real<lower=0,upper=1> wI[Nsub];                 '
   '    real aI[Nsub];                                  '
   '    real bI[Nsub];                                  '
   '}                                                   '
  
   'model {                                                        '
   '    int n;                                                     '
   '    real Pr;                                                   '
   '    real pay1;                                                 '
   '    real pay2;                                                 '
   '    real w1;                                                   '
   '    real w2;                                                   '
   '    mw1 ~ uniform(1,100);                                      ' 
   '    mw2 ~ uniform(1,100);                                      '  
   '    ma ~ normal(1,10);                                         '
   '    mb ~ normal(0,10);                                         '
   '    stda ~ uniform(.01,10);                                    ' 
   '    stdb ~ uniform(.01,10);                                    '
   '    wI ~ beta(mw1,mw2);                                        ' 
   '    aI ~ normal(ma,stda);                                      '
   '    bI ~ normal(mb,stdb);                                      '
   '    for (j in 1:Nsub) {                                        '
   '        for (k in 1:Nt)    {                                   '
   '            n <- (j-1)*Nt + k;                                 '
   '            w1 <-     wI[j]*(1-PD[n,3])+PD[n,3]*PD[n,2];       '
   '            w2 <- (1-wI[j])*(1-PD[n,3])+PD[n,3]*(1-PD[n,2]);   '
   '            pay1 <- w1*PM[G[n],1]+w2*PM[G[n],3];               '
   '            pay2 <- w1*PM[G[n],5]+w2*PM[G[n],7];               '
   '            Pr <- aI[j]*((pay1-pay2)) + bI[j]*PD[n,3];         '
   '            PD[n,1] ~ bernoulli_logit(Pr) ;                    '
   '        }                                                      '  
   '    }                                                          '
   '}                                                              '
};


  fit = stan('model_code',MervModel,'data',LData,'chains',4,'iter',10000,'warmup',5000,'thin',1);
% fit2 = stan('fit',fit1,'data',LData,'iter',10000,'chains',4);


 print(fit)
% fit.traceplot  (don't use!)

% return a struct with all parameters when none specifically requested
% samples = fit.extract('permuted',false);  
% AM = samples.ma;    % 16000 x 6 matrix: 4*4000 = 16000 samples of 6 alpha parms 
% BM = samples.mb;     % 16000 x 6 matrix: 4*4000 = 16000 samples of 6 beta parms
% mean(AM)
% mean(BM)
% log_lik = fit.extract.log_lik;
% mstan.waic(log_lik)
% w1 = samples(1).w; 
% b2 = samples(2).b;
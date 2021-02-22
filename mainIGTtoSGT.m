% main program used to estimate parameters
% generalize from IGT to SGT
clear
clc

load IGT4001to4002.txt   % Iowa gambling task
load SGT4001to4002.txt   % Suzhou gambling task

% 240 rows by 5 columns
% first 120 rows is subj1
% second 120 rows is sub2

% the data set has 120 obs from 2 subjects 2*120 = 240
% to pick up first subject use rows 1 - 120
% to pick up second subject use rows 121-240

% col=1 is originally trial no., which we don't need
% col=2 is deck
% col=3 is reward
% col=4 is loss
% col=5 is subj #

% now we eliminate col=1, and only pick 2,3,4


nt = 120;
sub = 1;

if sub == 1
    dataIGT = IGT4001to4002(1:nt,[2 3 4]);   %  first of two subjects
    dataSGT = SGT4001to4002(1:nt,[2 3 4]);   %  same first of two subjects
else
    % to pick up second subject use
    dataIGT = IGT4001to4002(nt+(1:nt),[2 3 4]);   %  second of two subje
    dataSGT = SGT4001to4002(nt+(1:nt),[2 3 4]);   % same second subj
end

dataE = dataSGT;   % Est 
dataG = dataIGT;   % Gen test


% pick the columns needed for model 
% first col = deck chosen
% sec col = win amt
% third col = lose amt

%  IGT deck number interpretation (in notes)
% 1 = deck  B, 100, 250
% 2 = deck A,  100, 1250
% 3 = deck C,  50 ,50
% 4 = deck D   50, 250

%  SGT (using notes)
% 1 = deck  A, 100, 525
% 2 = deck B,  50, 325
% 3 = deck C,  525 , 100
% 4 = deck D   325, 50

options = optimset('fminsearch');
% options =optimset(options,'Display','iter');

parm0 = [0  0  0];
% parm0 = randn(3,1);
% Q0 = zeros(4,1);

[parmE, ChiE1] = fminsearch(@(parm) FitEVL1(parm,dataE) , parm0,options);

'ChiE1'
ChiE1

y = dataE(1:nt,1);
parm0 = [mean(y==1) mean(y==2) mean(y==3) mean(y==4)];
'G^2  binomial'
bin1(parm0,dataE)
[parmB, ChiB1] = fminsearch(@(parm) bin1(parm,dataE), parm0, options);


% convert matlab real values into model bounded values

eta = exp(parmE(1));   % learning rate
eta = eta./(1+eta);

w = exp(parmE(2));    % weight for loss
w = w./(1+w);

c = exp(parmE(3));    % exploration/exploitation

disp('EVL parm')
disp('  learning   losswgt    choice')
disp([eta w c])

% 'Baseline Chi Sq minus EVL Chi Sq '
disp('Chi improvment Base Chi - EVL Chi fit')
disp(' positive values favor EVL over baseline')
disp(ChiB1-ChiE1)

% generalization to SGT
parmR = [1 1 1 1]/4;
ChiB2 = bin1(parmR,dataG);
'Generalization Chi'
ChiE2 = FitEVL1(parmE,dataG)

% 'Base Chi Sq minus EVL Chi Sq '
disp('Base Chi - EVL Chi generalization')
disp(' positive values favor EVL over baseline')
disp(ChiB2-ChiE2)



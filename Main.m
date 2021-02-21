% Merv Belief based best response model
clc
clear

load Merv1
% PD  88 subj, 40 trials per subj
% 40 trials, 28 no bonus , 12 bonus
% N = 40*88 = 3520, 28*88 = 2464
% N x 10 data matrix
%    1    2    3       4      5        6        7      8     9   10
%  subj trial matrix choice opponent oppchoice bonus payoff RT   Gnumb
ns = 88; % no subj for individual fits

% [By c] = find(PD(:,7)==2);  %  bonus trials
% [Bn c] = find(PD(:,7)==1);  %  not bonus
[G1, c ] = find(PD(:,10)==1); %  PD only   18 trials, 12 bonus 6 no bonus
PD = PD(G1,:);


PD(:,4) = 2 - PD(:,4);   % 1 = defect, 0 = coop
N = size(PD,1);  % 18*88 = 1584 only PD, 28*88 = 2464 no bonus,  or 40*88 = 3520 all
nt = N/ns;  % no trials per subj  18 PD only or 28 (no bonus) or 40 (all)

load Matrices   % payoff matrices used in games
% PMall 6 matrices x 8 cols x 5 games
% 1=PD,2=CH,3=IS,4=BS,5=FR : types of games
% cols
% M Y  M Y  M Y  M Y
% D D  D C  C D  C C
% 1 1  1 2  2 1  2 2  num code for D,C
% col's contain payoffs for a game

np = 3;   % no parms: wgt per game (6), choice parm (1), bonus effect (1)
ParmV = zeros(ns,np);
ChiV = zeros(ns,1);

nr = 5;  % no of random start parms
lb = zeros(np,1);
% ub = [ones(5,1); 10.*ones(np-5,1)];
ub = [1 ; 10*ones(np-1,1) ]';   % for PD model only

if 1 == 2
    
    for sub = 1:ns
        sub
        RepParm = zeros(nr,np);
        RepChi = zeros(nr,1);
        for  rep = 1:nr
            % Parm0 = randn(np,1);
            Parm0 = rand(np,1);
            %  Parm0 = [.4 .2 10]';
            indx = (sub-1).*nt + (1:nt);  % individual subject
            PS = PD(indx,:);
            options = optimset('MaxFunEvals',10000,'MaxIter',10000,...
                'TolX',1e-10);
            [param, Chi] = fmincon(@(Parm) FitMerv3(Parm, PS, PMall, nt),...
                Parm0, [],[],[],[],lb,ub,[],options);
            [rep Chi/1000]
            RepParm(rep,:) = param;
            RepChi(rep) = Chi;
        end % rep
        [ mv mindx ] = min(RepChi);
        ParmV(sub,:) = RepParm(mindx,:);
        ChiV(sub) = mv;
    end % subj
    
     %  save MFit3 ParmV ChiV
    
else
    
    load MFit3
    Comp = zeros(ns,1);
    Pred = zeros(ns,6,3);
    % six games, 2 bonus trials per game,
    % 1 no bonus trial per game
    
    for sub = 1:ns
        indx = (sub-1).*nt + (1:nt);
        PS = PD(indx,:);
        Parm = ParmV(sub,:);
        
        [Pr, Chi] = FitMerv3p(Parm,PS,PMall,nt);
        Pred(sub,:,:) = Pr;
        
    end
    
    M = squeeze(mean(Pred,1));
    disp(M)  % for avg of indiv fit
    mean(M)
    
end  % if statement


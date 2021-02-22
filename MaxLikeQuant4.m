% Max Like Est quant Theory
% four parameter with gamma
clc
clear

% nsub = 100
% 33 choices per subj = 16 gambles played twice plus ne gamble at start
% data is a nsub*33 X 12: g l trialpay stagepay  
%    wptt wptn wpnt wpnn lptt lptn lpnt lpnn
% avg log like ranges from -104 to -37

load BarkanData  % contains a file called: data
size(data)
data = reshape(data,100,33,12);
load MaxLikeFitsQ4
  
  nsub= 100;   % 100
  SBF = []; 
  
  
% options = optimset('fminsearch'); 
% oldopts = optimset('fminsearch'); 
% options = optimset(oldopts,'Display','iter');
 options = optimset('MaxFunEvals',1000,'MaxIter',1000,'TolX',1e-5);
 
w = 20;                    % 10  even number
W = 2*w+1;                 % 41 values
W = ((1:W)'-(w+1))/w;      % -1 to 0   to  +1  inc .05
UtP = .9 + .5*W;           % .4 to 0.9 to 1.4  inc .025  
LP = (1.5+W);              % .5 to 1.5 to 2.5  inc .05
MemP = .5*(1+W);           %  0 to 0.5 to 1.0  inc .025
GamQ = 5*W;                % -5 to 0.0 to  5   inc .25    
GamP = 2*(1+W);            %  0 to  2  to  4   inc .10 
Parms = [UtP LP MemP GamQ GamP];

 
% [v1,v2,v3,v4]=ind2sub([41 41 41 41 41],ML(sub,2))
% LQ = Lquant(Sdata,Parms(v1,1),Parms(v2,2),Parms(v3,3),Parms(v4,4))
 
 
SML = [];
  mp = []  ;
if 1==1
LQm = [];
for sub= 1:nsub
    sub
    Sdata=squeeze(data(sub,:,:));
    [v1,v2,v3,v4]=ind2sub([41 41 41 41 41],ML(sub,2));
    parm0 = [Parms(v1,1),Parms(v2,2),Parms(v3,3),Parms(v4,4)];
    LQ = Lquant(parm0,Sdata);
    LQm = [LQm ; LQ];
    [parm,MLs] = fminsearch(@(parm) Lquant(parm, Sdata), parm0, options);
    [MLs LQ]
    mp = [mp ; parm];
   
    SML = [SML ; [sub MLs]];
   % save MLfileQ4 SML
end

mean(SML)
median(SML)
std(SML)
mean(LQm)

end



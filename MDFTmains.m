% June 21 2017
% main program
% estimate parms and generate predictions for MDFT
%  fits the choice probabilities for n alternative choice
% you can use m attributes but distfunct only works with 2
% uses the subject controlled, self terminating stopping rule
% uses distance function described in Hotaling article


clear
clc
Ns = 10000;  % no. of simulations per condition

sw =2;    % 1 = estimate parameters , 2 = skip fit , generate predictions

% describe the choice options
% three options on two attributes
% 3 x 2 x N matrix  , N = no. of choice sets
% 
%             op A    opt B    opt C
MM = zeros(3,2,10);
MM(:,:,1) = [ 1   3;  3  1;  .9  3.1 ];  % sim
MM(:,:,2) = [ 1   3;  3  1; 1.1  2.9 ];  % sim
MM(:,:,3) = [ 1   3;  3  1; 1.75 2.25];  % com
MM(:,:,4) = [ 1   3;  3  1; 2    2   ];   % comp
MM(:,:,5) = [ 1   3;  3  1; 2.25 1.75] ;  % comp
MM(:,:,6) = [ 1   3;  3  1;  .5  2.5 ];   % attraction
MM(:,:,7) = [ 1   3;  3  1;  1.1 2.5 ];   % attraction
MM(:,:,8) = [.5  .5; .7 .7;   2    2 ];  % dom
MM(:,:,9) = [.75 .75; 1  1; 1.25 1.25];  % dom
MM(:,:,10) =[.5  .5;  1  1; 2    2   ];  % dom


% hypothetical results for 10 these choice sets
% do not take too seriously
% D is a matrix of choice proportions
% each row is a choice set, each col. is an alternative

D = [ .3   .4   .3 ;
    .3   .4   .3 ;
    .3  .34  .36;
    .3   .3   .4 ;
    .34  .3  .36;
    .8   .2    0 ;
    .8   .2    0 ;
    0    0    1 ;
    0    0    1 ;
    0    0    1];

N = 100;  % no. observations per choice set 
% alternatively enter a new matrix called FD containing the choice frequencies
FD = N*D; 

na = size(D,2);   % no of options in choice set

C3 = -1/(na-1)*ones(na,na);   % C matrix for na options
C3 = C3 - diag(diag(C3)) + eye(na);


%      dis wgt  dist phi1 dist ph2  sig    thresh    w=attention to attribue 1
% x0 = [log(12) log(.022) log(.05) log(1) log(17.5) .5];
 % x0 = [6.2197    0.0098    0.0482    1.0007   21.1196 .50]; % SSE = .1176
% x0 = [6.2197    0.00      0.0482    1.0007   21.1196 .50];
% x0 = log(x0);    % real valued matlab parms get exponentiated to positive values
  x0=[  1.8540   -4.5461   -3.0421    0.0007    3.1017   -0.7082]; % .0953

if sw==1  % fit data
    
    options = optimset('Display','iter','Maxiter',100);
    
    [x sse] = fminsearch(@(x) fitMDFTs(x,D,FD,MM,C3,Ns), x0,options);
    disp('matlab parms')
    disp(x)
    disp('sse')
    disp(sse)
    wgt = exp(x(1))
    phi1 = exp(x(2))
    phi2 = exp(x(3))
    sig2 = exp(x(4))
    theta1 = exp((5))
    
else  % by pass fitting
    [sse, P3, TV] = fitMDFTs(x0,D,FD,MM,C3,Ns);
    
    disp('sse')
    disp(sse)
    disp('         target data           model predictions')
    disp([D P3 TV])
    
end

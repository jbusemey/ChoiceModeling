function [P3, T] = simMDF(G3,C3,M3,w,theta1,sig2,Ns)
% function [P3 T] = simMDF(G3,C3,M3,w,theta1,sig2,Ns)

% S3 is the n x n gamma feedback matrix
% C3 is a contrast matrix 
% M3 is a n x m matrix of values of each option (row) on each attribute (col)

h = 1; hh = sqrt(h);     % time unit (can be changed)
na = size(w,1);
Na = 1:na; Na = Na' ;
% don't need to change
nd = size(M3,1);         % number of options
V3 = C3*M3*w;            % mean valence
P3 = zeros(1,nd);        % initial value of counter
T = 0;                   % initial value for time
 
% multipd = makedist('Multinomial','Probabilities',w) ;
Ind = 0;
for ns = 1:Ns
    B = 0; t=0; P = zeros(nd,1);   
    while (B < theta1)
     % WV = random(multipd,1,1);
        W = w(1) > rand; 
        W = [W ; (1-W)];    % pick an attribute
     %  W = ( WV == Na);
       E3 = C3*M3*(W-w.*h) + sig2*C3*randn(nd,1);  % compute noise 
       P = P + G3*P*h + V3*h + E3*hh;   % accumulate
       [B,Ind] = max(P);      % find max
       t = t+h;               % track time
  %     [n t]
    end % one simulation
    P3 = P3 +  [Ind ==1  Ind==2 Ind==3];
    T = T + t;
end % all N  simulations

P3 = P3./Ns;
T = T./Ns;

% if nd ==3
% else
%     P3 = [sum(P3 ==1) ; sum(P3==2); sum(P3==3); sum(P3==4)]/N;
% end


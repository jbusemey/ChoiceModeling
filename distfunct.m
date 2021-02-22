
function [D, EG] = distfunct(X,wgt,phi1,phi2)
% D = distfunct(X,wgt,phi1,phi2);
% distance functiton for lateral inhibition
% only for 2 - attributes


T = [-1 1 ; 1 1]./sqrt(2);
n = size(X,1);
W = diag([1 wgt]);  % wgt second dim more to emphasize dominance
D = zeros(n,n);
for i = 1:n
    for j = 1:n
        DV = (X(i,:)-X(j,:))';
        DV = T*DV;
        D(i,j) = DV'*W*DV;
    end
end
 D = - phi2*exp(-phi1*D.^2);      % gamma
% D = eye(n) - phi2*exp(-phi1*D);

EG = eig(D) ;
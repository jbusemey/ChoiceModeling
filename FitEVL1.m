function chi = FitEVL1(parm,data)
% function chi = FitEVL1(parm,data)

chi= 0 ;	    %initialize

y = data(:,1);
x = data(:,2:3);
nt = size(y,1);

eta = exp(parm(1));
eta = eta./(1+eta);

w = exp(parm(2));
w = w./(1+w);

c = exp(parm(3));

Q = zeros(4,1);

%-----------------------------------MODEL-----------------------
%begin of training
for tt = 2:nt
    t = tt-1;
    yt = y(t,:) == [1 2 3 4]';   % current choice just made this choice
    ytt= y(tt,:)== [1 2 3 4]';  % the next choice one step ahead we wish to predict
    xt = x(t,:);    % current payoff just experienced; one row in x, for trial t
    
    win  = abs(xt(1))/100;                                          % Eldad
    loss = abs(xt(2))/100;                                          % determine loss and rescale it
    
    Q = Q + eta * yt.* ((1-w)*win - w*loss - Q);
     th = (t/10)^c;    % sensitivity
    
    % th = c;
    
    s = exp(th*Q)+ 0.0000000001;  % strength
    p = s./sum(s);    % prob of choice for each deck
    
    pp = .0001 + .9998*p;
    chi = chi + log(pp)'*ytt;
    
end
%end of training
%-----------------------------------------MODEL----------------------------------
chi = -2*chi;  % -2 times log likelihood


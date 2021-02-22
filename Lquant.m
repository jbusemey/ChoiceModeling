function L=Lquant(parm,Sdata)
% function L=Lquant(parm,Sdata)

t = 1;
% ut = .4 + (1/(1+exp(-parm(1))));
% Lp = .5 + 2*(1/(1+exp(-parm(2))));
% mem=  (1/(1+exp(-parm(3))));
% gam= parm(4);   
ut = parm(1);
Lp = parm(2);
mem = parm(3);
gam = parm(4);

SU = [1 1 1 1]'/2;
SW = [1 1 0 0]'/sqrt(2);
SL = [0 0 1 1]'/sqrt(2);

M = [1 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 0];

L = 0;

% PM = zeros(33,9);

for k=1:33
    
    Wt = sum(Sdata(k,5:8));
    Lt = sum(Sdata(k,9:12));
    S = Wt.*SW + Lt.*(SL);
    
    x = Sdata(k,1:2)/100;
  
    %  PM(k,1) = (x(1)+x(2))/2;
    
    euw1 = (2*x(1)).^ut;
    euw2 = x(1)+x(2);
    euw2 = (euw2>=0)*euw2.^ut -(euw2<0)*Lp*abs(euw2).^ut;
    euw = .5*(euw1+euw2)-x(1).^ut;
    
    eul1 = x(2)+x(1);
    eul1 = (eul1>=0)*eul1.^ut -(eul1<0)*Lp*abs(eul1).^ut;
    eul2 = -Lp*abs(2*x(2)).^ut;
    eul = .5*(eul1+eul2)-(-Lp*abs(x(2)).^ut);
    
    mu1 = 2./(1+exp(-euw))-1;
    mu2 = 2./(1+exp(-eul))-1;
       
    Ha1 = (1/sqrt(1+(mu1.^2)))*kron([1 0 ; 0 0],[mu1 1; 1 -mu1]);
    Ha2 = (1/sqrt(1+(mu2.^2)))*kron([0 0 ; 0 1],[mu2 1; 1 -mu2]);
    Ha = Ha1 + Ha2;
    
    Hc = (-gam/sqrt(2))*( kron([1 1;1 -1],[1 0; 0 0])...
        + kron([-1 1;1 1],[0 0; 0 1]) );
    
    H = Ha + Hc;
    U = expm(-1i*t*(pi/2)*H);
    
    Sp = U*SU;
    Pp = Sp'*M*Sp;
    
    Sf = U*S;
    Pf = Sf'*M*Sf;
    
    Ptt = mem*Pp + (1-mem)*Pp*Pf;
    Ptn = (1-mem)*Pp*(1-Pf);
    Pnt = (1-mem)*(1-Pp)*Pf;
    Pnn = mem*(1-Pp) + (1-mem)*(1-Pp)*(1-Pf);
    PP = [Ptt Ptn Pnt Pnn]';

    % PM(k,2:9) = [Wt*PP' Lt*PP'];
    
    eps = .00001;
    PP = eps + (1-eps)*PP;
    lnP = log(PP);
    D = (Wt*Sdata(k,5:8)+Lt*Sdata(k,9:12));
    L = L + D*lnP;
        
end
L = -2*L;

% PM

% Analyze results
%
clear
clc


 
 load samplesI3_6;
 W = ParmsH(:,1)./(ParmsH(:,1)+ParmsH(:,2));
 ParmsH = [W ParmsH(:,3:6)];
 V = ParmsH;
%  W = ParmsH(:,1)./(ParmsH(:,1)+ParmsH(:,2));
%  ParmsH = [W ParmsH(:,3:6)];
%  ParmsH = reshape(ParmsH,4,10000,5);
%  ParmsH = squeeze(ParmsH(4,:,:));
%  V =  ParmsH;
%  

fig = 3;

% Use to compute R
if fig == 6
    nr = 4;
    ns = size(V,1)/nr;
    np = size(V,2);
    V=reshape(V,ns,nr,np);
    V = permute(V,[2 1 3]);
else
    ns = size(V,1);
    np = size(V,2);
end


k = 0;  


if 1==fig
    subplot(3,1,1)
    plot(V(:,k+1))
    ylabel('w slope');
    subplot(3,1,2)
    plot(V(:,k+2))
    ylabel('a int');
    xlabel('Sample Number');
    subplot(3,1,3)
    plot(V(:,k+3))
    ylabel('b slope');
    
    
elseif 2==fig
    movavg = filter( (1/50)*ones(50,1), 1, V);
    
    subplot(3,1,1)
    plot(movavg(:,1))
    xlabel('Number of samples')
    ylabel('Means of the intercept');
    subplot(3,1,2)
    plot(movavg(:,2))
    xlabel('Number of samples')
    ylabel('Means of the slope');
    subplot(3,1,3)
    plot(movavg(:,3))
    xlabel('Number of samples')
    ylabel('Means of the intercept');
    
    
elseif 3==fig
    
    npt=100;
    pts = (0:npt);
    pp = zeros(1,npt+1);
    ci = 6:96;  % 90%
    
    subplot(3,1,1)
    xi = prctile(V(:,k+1),pts);
    [f, xi] = ksdensity(V(:,k+1),xi);
    fm = max(f);
    px = pp; px(ci)=.05*fm; px(51) = fm;
    plot(xi,f,'-',xi,px,'-')
    xlabel('w');
    ylabel('Posterior Probability')
    
    subplot(3,1,2)
    xi = prctile(V(:,k+2),pts);
    [f, xi] = ksdensity(V(:,k+2),xi);
    fm = max(f);
    px = pp; px(ci)=.05*fm; px(51) = fm;
    plot(xi,f,'-',xi,px,'-')
    xlabel('a');
    ylabel('Posterior Probability')

    subplot(3,1,3)
    xi = prctile(V(:,k+3),pts);
    [f, xi] = ksdensity(V(:,k+3),xi);
    fm = max(f);
    px = pp; px(ci)=.05*fm; px(51) = fm;
    plot(xi,f,'-',xi,px,'-')
    xlabel('b');
    ylabel('Posterior Probability')
    
    disp('posterior means')
    disp('Wgt a b stda stdb')
    disp(mean(V))
    
elseif 4==fig
    
    subplot(3,1,1)
    hist(V(:,k+1))
    xlabel('Cat slope');
    subplot(2,3,2)
    hist(V(:,k+2))
    xlabel('C&D0 int');
    
    subplot(3,1,2)
    hist(V(:,k+3))
    xlabel('C&D0 slope');
    subplot(2,3,4)
    hist(V(:,k+4))
    xlabel('D2 slope');
    
    subplot(3,1,3)
    hist(V(:,k+5))
    xlabel('C&D1 int');
    subplot(2,3,6)
    hist(V(:,k+6))
    xlabel('C&D1 slope');
    
elseif 5==fig
    timedata = iddata(V);
    % get(timedata)
    y1 = timedata.y(:,1);
    y1 = (y1-mean(y1))./std(y1);
    nt = size(y1,1);
    armodel = ar(y1,5)
    armodel.a
    armodel.NoiseVariance
    
elseif 6==fig
    M = squeeze(mean(V,2));
    MM = ones(nr,1)*mean(M);
    BD = M-MM;
    BD = ns*diag(BD'*BD)/(nr-1);
    
    WD = zeros(size(V,3),1);
    for nreps = 1:nr
        D = squeeze(V(nreps,:,:))-ones(ns,1)*M(nreps,:);
        WD = WD + diag(D'*D);
    end
    WD = WD/(nr*(ns-1));
    
    SD = ((ns-1)*WD + BD)/ns;
    R = sqrt(SD./WD)
    
elseif 7 == fig
    % Exp 1
    %   1       2      3      4       5     6
    %   cat   C&D0 1 C&D0 2   D2   C&D1 1  C&D1 2  intercept
    %   7       8      9      10     11     12
    %   cat   C&D0 1 C&D0 2   D2   C&D1 1  C&D1 2   slope
    %  Exp 2
    %   1       2      3      4       5
    %   cat   C&D0 1 C&D0 2   D2     D3    intercept
    %   6       7      8      9      10
    %   cat   C&D0 1 C&D0 2   D2     D3    slope
    
    var1 = 4; var2 = 5;
    var1 = [ var1 var1+5];
    var2 = [ var2 var2+5];
    VD = V(:,var1)-V(:,var2);
    
    if 1 == 2
        vv = 1;  % 1 = int , 2 = slope
        np=100;
        pts = (0:np);
        % ci = 6:96;  % 90%
        % ci = 3:99;  % 96%
        ci = 4:98;    % 94%
        pp = zeros(1,np+1);
        xi = prctile(VD(:,vv),pts);
        dx = [diff(xi) 0];
        [f, xi] = ksdensity(VD(:,vv),xi);
        fm = max(f);
        px = pp; px(ci)=.01*fm; px(51) = fm;
        plot(xi,f,'-',xi,px,'-')
        xlabel('Difference');
        ylabel('Relative Frequency')
        title('C&D0 1 minus C&D1 1 Int')
        
    else
        plot(VD(:,1),VD(:,2),'.',0,0,'o',mean(VD(:,1)),mean(VD(:,2)),'*')
        xlabel('Intercept')
        ylabel('Slope')
    end
    
end

wI = ParmsI(:,1:88);
aI = ParmsI(:,89:(2*88));
bI = ParmsI(:,((2*88)+1):3*88);


% reshape(mean(squeeze(mean(simdat,2))),np/2,2)

% std(V)
% corr(V)

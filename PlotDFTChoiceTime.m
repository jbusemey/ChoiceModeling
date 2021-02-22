clear
clc
TV = (0:.01:3)';
parms = [1 1 0 -.1];
[pr, M, pdf, cdf] = DFT2ALTCRT(parms, TV);


plot(TV,pdf(:,1),'-',TV,pdf(:,2),'-')
title('Distribution of Choice Time')
legend('A','B')
xlabel('Choice Time')
ylabel('Probability density')



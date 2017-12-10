% To check whether likelihood is decreasing 
%

clear all
iterN = 28;

rLL1 = zeros(iterN,1);
rLL2 = zeros(iterN,1);
rLL = zeros(iterN,1);
rLLdiff = zeros(iterN,1);

for iterI = 1 : iterN
    load(sprintf('0_iter%d.mat', iterI))
    rLL1(iterI) = tempLL1;
    rLL2(iterI) = tempLL2;
    rLL(iterI) = loglik;
    rLLdiff(iterI) = comdiff;
end

subplot(4,1,1); plot(rLL(5:end))
subplot(4,1,2); plot(rLLdiff(5:end))
subplot(4,1,3); plot(rLL1(5:end))
subplot(4,1,4); plot(rLL2(5:end))

%%
clear all
% iterN = [59 48 25 48 51 53 36 39 38 56 36 36 46 34 62 48 55 53 36 45];
iterN = [55 52 54 47 38 35 48]
% iterN = [24 25 28 44 29  23 37 57 22 27   26 22 48 24 19   27 21 35];
maxIterN = max(iterN);
runN = 7;


rLL1 = NaN(maxIterN,runN);
rLL2 = NaN(maxIterN,runN);
rLL = NaN(maxIterN,runN);
rLLdiff = NaN(maxIterN,runN);
lastrLL = zeros(runN,1);

totalthetaE1 = zeros(4,4,runN);
totalthetaE2 = zeros(4,4,runN);
totalsigmaE1 = zeros(4,4,runN);
totalsigmaE2 = zeros(4,4,runN);

for runNi = 1 : runN
    for iterI = 1 : iterN(runNi)
        load(sprintf('%d_iter%d.mat', runNi, iterI))
        rLL1(iterI,i) = tempLL1;
        rLL2(iterI,i) = tempLL2;
        rLL(iterI,runNi) = loglik;
        rLLdiff(iterI,runNi) = comdiff;
    end
    lastrLL(runNi) = loglik;
    totalthetaE1(:,:,runNi) = thetaE1;
    totalthetaE2(:,:,runNi) = thetaE2;
    totalsigmaE1(:,:,runNi) = sigmaE1;
    totalsigmaE2(:,:,runNi) = sigmaE2;
end

startP = 1;
% subplot(4,1,1); plot(rLL(startP:end,:))
% subplot(4,1,2); plot(rLLdiff(startP:end,:))
% subplot(4,1,3); plot(rLL1(:,startP:end))
% subplot(4,1,4); plot(rLL2(:,startP:end))
% 
subplot(2,1,1); plot(rLL(startP:end,:))
subplot(2,1,2); plot(rLLdiff(startP:end,:))

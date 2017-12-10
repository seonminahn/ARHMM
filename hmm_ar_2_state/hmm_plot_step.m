% function [e, a, max_lk, c1, c2, st] = hmm_plot_step(data)
% this runs the B-W routine 20 times and returns the max likelihood values
% it plots the data and backsamples to get teh state probabilities at each
% position

% this is where Luan's initial values are loaded
% % % load('cov3.mat');

function [sa, a, thetaE1, thetaE2, sigmaE1, sigmaE2, means, max_lk, dataFilled] = hmm_plot_step(data, run)
% clear all
% load('4DdataWtime_norm.mat');   % fourvariablesCEL
% data = fourvariablesCEL';
% data = fliplr(data);
% data(1,:) = 2*abs(data(1,:) - data(1,1));


% figure;
% for i = 2:size(data,1)
%     subplot(4, 1, i-1)
%     plot(data(1, :), data(i, :))
%     ylim([min(data(i,:)), max(data(i,:))]);
%     title('Data')
% end

% run B-W a bunch of times and get teh max log likelihood and parameters
max_lk = -realmax;
iter = 1;
[nvars, ~] = size(data);
nvars = nvars - 1;
rand('seed', sum(clock)+run);
while iter <= 1
    % initailize with Luan's results
    % it works with other reasonable starting values also

%     load('parameterE.mat'); % sigmaE1, sigmaE2, thetaE1, thetaE2
%     ithetaE1 = thetaE1;
%     ithetaE2 = thetaE2;
%     isigmaE1 = sigmaE1;
%     isigmaE2 = sigmaE2;

    % Start with random values
    ithetaE1 = genE(nvars);
    ithetaE2 = genE(nvars);
    ithetaE1 = ithetaE1/norm(ithetaE1,1);
    ithetaE2 = ithetaE2/norm(ithetaE2,1);
%     tempE = genE(dataDim); isigmaE1 = tempE*tempE';
%     tempE = genE(dataDim); isigmaE2 = tempE*tempE';
    tempE = rand(nvars,nvars) + 4*eye(nvars) - 0.5; isigmaE1 = tempE*tempE'/100;
    tempE = rand(nvars,nvars) + 4*eye(nvars) - 0.5; isigmaE2 = tempE*tempE'/200;

    [loglik, tsa, ta, tthetaE1, tthetaE2, tsigmaE1, tsigmaE2, tmeans, dataFilled]=n_baum_welch(iter, data, ithetaE1, ithetaE2, isigmaE1, isigmaE2);
    if loglik > -realmax
        disp(['iteration ' num2str(iter)])
        if loglik > max_lk
            max_i = iter;
            max_lk = loglik;
            % % %             a = a1;
            % % %             e = e1;
            % % %             c1 = cv1;
            % % %             c2 = cv2;
            sa = tsa;
            a = ta;
            thetaE1 = tthetaE1;
            thetaE2 = tthetaE2;
            sigmaE1 = tsigmaE1;
            sigmaE2 = tsigmaE2;
            means = tmeans;
        end
        iter = iter + 1;
    end
end

fprintf('\n');
disp(['Log likelihood = ', num2str(max_lk), ' iteration = ', num2str(max_i)]);
disp('Initial Matrix:');
disp(sa);
disp('Transition Matrix:');
disp(a);
disp('ThetaE1:');
disp(thetaE1);
disp(eig(thetaE1));
disp('ThetaE2:');
disp(thetaE2);
disp(eig(thetaE2));
disp('SigmaE1:');
disp(sigmaE1);
disp('SigmaE2:');
disp(sigmaE2);
disp('constant');
disp(means)

% % sample to get state probs
% % this is the slow part
% nSample = 1000;
% st=backsample(data, nSample, sa, a, thetaE1, thetaE2, sigmaE1, sigmaE2, means);
% 
% % % plot the sampled state probs twice
% % figure;
% % for i = 1:2
% %     subplot(2,1,i)
% %     plot(sum(st==i)/nSample)
% %     title(['Probability of Being in State ' num2str(i)]);
% %     ylabel('probability');
% % end
% % 
% % figure;
% % for i = 1:2
% %     subplot(2,1,i)
% %     plot(sum(st==i)/nSample, 'o')
% %     title(['Probability of Being in State ' num2str(i)]);
% %     ylabel('probability');
% % end
% 
% % figure;
% % for i = 1:2
% %     subplot(2,1,i)
% %     stSum = sum(st==i)/1000;
% %     depthI = fliplr(x(1,2:end)/2);
% %     plot(depthI, stSum(2:end), 'o')
% %     title(['Probability of Being in State ' num2str(i)]);
% %     xlabel('depth')
% %     ylabel('probability');
% % end
% 
% % figure;
% % for i = 1:2
% %     subplot(2,1,i)
% %     stSum = sum(st==i)/1000;
% %     plot(data(1,:)/2, stSum, 'o')
% %     title(['Probability of Being in State ' num2str(i)]);
% %     xlabel('Time')
% %     ylabel('probability');
% % end
% 
% save('hmm.mat');



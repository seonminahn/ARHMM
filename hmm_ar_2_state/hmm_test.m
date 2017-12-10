% this just a wrapper that rund the B-W routine and save the data
% it runs 5 times and makes zip files with the plots and output for later
% examination. This way you can go get a snack while it runs. :-)

function hmm_test(sampleName, sampleLoc, resultsLoc, runFor, figurePlot)

tic
% clear;
max_log = -realmax;
max_run = 0;
for run = runFor
%     close all
%     zip_file = ['fit_', num2str(run), '.match.zip'];
    
    % Original data set
    % % %     % data is in d2.mat
    % % %     % Luan's results in cov3.mat
    % % %     load('d2.mat');
    % % %     load('cov3.mat');
    
    % Load the normalized data and flip it.
%     load('4DdataWtime_norm.mat');   % fourvariablesCEL
%     data = fourvariablesCEL';
%     data = fliplr(data);
%     data(1,:) = 2*abs(data(1,:) - data(1,1));
    
    % Load simulated data
    load(sampleLoc)
    data = sampleD';
    nvars = size(data,1)-1;
 
    %%%%% st is added by Seonmin
    % % %     [e2, a2, max_lk2, c1, c2, st] = hmm_plot_step(data);
    [sa, a, thetaE1, thetaE2, sigmaE1, sigmaE2, means, max_lk2, dataFilled] = hmm_plot_step(data, run);
    
    if max_lk2 > max_log
        max_log = max_lk2;
        max_run = run;
%         fid2 = fopen('max_run.txt', 'w');
%         fprintf(fid2, 'max log likelihood = %f run %d\n', max_log, max_run);
%         fclose(fid2);
    end
    
%     print(1, '-dpdf', 'data.pdf');
%     print(2, '-dpdf', '2_state_fit.pdf');
%     print(3, '-dpdf', '2_state_fit_points.pdf');

    fprintfFormat = '';
    for i = 1 : nvars
        fprintfFormat = sprintf('%s%%8.4f ', fprintfFormat);
    end
    fprintfFormat = sprintf('%s\n', fprintfFormat);

    fid = fopen(sprintf('%s%s_fit_run%d-%d_results.txt', resultsLoc, sampleName, runFor(1), runFor(end)), 'a+');
    fprintf(fid, '\nrun %d\n', run);

    fprintf(fid, '\n2 State Model\n');
    fprintf(fid, '\n2 state log likelihood = %f\n', max_lk2);

    fprintf(fid, '\nInitial matrix:\n');
    fprintf(fid, '%f %f\n', sa);    
    
    fprintf(fid, '\nTransition matrix:\n');
    fprintf(fid, '%f %f\n', a');
    
    fprintf(fid, '\nPredicted Means\n');
    fprintf(fid, fprintfFormat, means);
%     fprintf(fid, 'Norm diff means: %f\n', norm(e2 - [m1; m2]));     

    cor = zeros(nvars, nvars, 2);
    cor(:,:,1)=thetaE1;
    cor(:,:,2)=thetaE2;    
    fprintf(fid, '\nPredicted Autocorrelation Matrix for State %d\n', 1);
    fprintf(fid, fprintfFormat, cor(:, :, 1)');    
    fprintf(fid, 'Eigenvalues %d\n', 1);
    fprintf(fid, fprintfFormat, real(eig(cor(:, :, 1))));
    fprintf(fid, fprintfFormat, imag(eig(cor(:, :, 1))));
    
    fprintf(fid, '\nPredicted Autocorrelation Matrix for State %d\n', 2);
    fprintf(fid, fprintfFormat, cor(:, :, 2)');
    fprintf(fid, 'Eigenvalues %d\n', 2);
    fprintf(fid, fprintfFormat, real(eig(cor(:, :, 2))));
    fprintf(fid, fprintfFormat, imag(eig(cor(:, :, 2))));

    cv = zeros(nvars, nvars, 2);
    cv(:,:,1)=sigmaE1;
    cv(:,:,2)=sigmaE2;    
    fprintf(fid, '\nPredicted Covariance Matrix for State %d\n', 1);
    fprintf(fid, fprintfFormat, cv(:, :, 1));
%     fprintf(fid, 'Norm diff state %d: %f\n', 1, norm(cov1 - cv(:, :, 1)));
    fprintf(fid, '\nPredicted Covariance Matrix for State %d\n', 2);
    fprintf(fid, fprintfFormat, cv(:, :, 2));
%     fprintf(fid, 'Norm diff state %d: %f\n', 2, norm(cov2 - cv(:, :, 2)));
    fclose(fid);
    
    save(sprintf('%s%s_results_dat%d.mat', resultsLoc, sampleName, run));
%     zip(zip_file, {'2_state_fit.pdf', '2_state_fit_points.pdf', 'data.pdf', ...
%         'fit_results.txt', 'dat.mat'});
    toc
end


disp(['max log likelihood = ', num2str(max_log), ' run ', num2str(max_run)]);
fid = fopen(sprintf('%s%s_fit_run%d-%d_results.txt', resultsLoc, sampleName, runFor(1), runFor(end)), 'a+');
fprintf(fid, '\nmax log likelihood = %f run %d\n', max_log, max_run);
fclose(fid);


% figure;
% stSum = sum(st==1)/1000;
% for i = 1:4
%     subplot(4,1,i)
%     plot(data(1,:), stSum*6-3, 'bo', data(1,:), data(i+1,:), 'r')
% end

% sample to get state probs
% this is the slow part
clearvars -except resultsLoc sampleName max_run
close all
load(sprintf('%s%s_results_dat%d.mat', resultsLoc, sampleName, max_run))
nSample = 1000;
data = dataFilled;
st=backsample(data, nSample, sa, a, thetaE1, thetaE2, sigmaE1, sigmaE2, means);
st(:,1) = [];
toc

% % plot the sampled state probs twice
% figure;
% for i = 1:2
%     subplot(2,1,i)
%     plot(sum(st==i)/nSample)
%     title(['Probability of Being in State ' num2str(i)]);
%     ylabel('probability');
% end
% 
% 
% figure;
% for i = 1:2
%     subplot(2,1,i)
%     stSum = sum(st==i)/1000;
%     depthI = fliplr(x(1,2:end)/2);
%     plot(depthI, stSum(2:end), 'o')
%     title(['Probability of Being in State ' num2str(i)]);
%     xlabel('depth')
%     ylabel('probability');
% end
% 
% figure;
% for i = 1:2
%     subplot(2,1,i)
%     stSum = sum(st==i)/1000;
%     plot(data(1,:)/2, stSum, 'o')
%     title(['Probability of Being in State ' num2str(i)]);
%     xlabel('Time')
%     ylabel('probability');
% end
if figurePlot == true
    figure;
    for i = 1:2
        subplot(2,1,i)
        plot(sum(st==i)/nSample, 'o')
%         currentAxis = axis; currentAxis(2) = size(st,2); axis(currentAxis)
        title(['Probability of Being in State ' num2str(i)]);
        ylabel('probability');
    end
    export_fig(sprintf('%s%s_fit_run%d.pdf', resultsLoc, sampleName, max_run)) 
end
save(sprintf('%s%s_results_dat%d.mat', resultsLoc, sampleName, max_run));

% If we generated a sample data, compare results with true states
sampleDataLoc = ['../samples/' sampleName '_data.mat'];
if exist(sampleDataLoc, 'file') == 2
    fid = fopen(sprintf('%s%s_fit_results.txt', resultsLoc, sampleName), 'a+');
    load(sampleDataLoc)
    ind1 = sample_statefull(2:end) == 1;
    ind2 = (sum(st==1)/1000 > 0.5)';
    sameIND = find(ind1 ~= ind2);
    disp(['Number of coinciding states: ' num2str(length(sameIND)) ' or ' num2str(length(sample_statefull)-1-length(sameIND))])
    fprintf(fid, '\nNumber of coinciding states: %d\n', length(sameIND));
    fclose(fid);
end

toc


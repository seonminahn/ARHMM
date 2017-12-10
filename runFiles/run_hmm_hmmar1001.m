clear all
currentFolder = pwd;

cd ..

% List sample name in sampleName
sampleName = [{'new4D_normalized_flipped_doubleTime'}];
solveMethod = [{'hmm_ar_2'}];
% If there are multiple samples, specify i to use the ith sample
i = 1;
j = 1;

% Folder path for sample files and result files
sampleLoc = ['../newData/' sampleName{i} '.mat'];
resultsLoc = ['../results_' solveMethod{j} '_state/'];

% Set runFor = 1:5 to run 5 times
% Depending on randomly selected initial values, the algorithm may converge
% to local optimal values or fail to initialization steps.
runFor = 1;

cd([solveMethod{j} '_state'])
mkdir(resultsLoc)
diary(sprintf('%s%s_diaryRun%d-%d.txt', resultsLoc, sampleName{i}, runFor(1), runFor(end)));
tic
hmm_test(sampleName{i}, sampleLoc, resultsLoc, runFor, 0)
toc
diary off
cd(currentFolder)
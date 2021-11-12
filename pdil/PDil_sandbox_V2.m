%AT 11/10/21, beginning to look at PDil analysis

%the V2 here is to differtiate from previous version of sandbox I'd started
%on my other laptop. Made it V2 mostly to prevent confusion in the future
%and prevent weird git things if I get the older version onto same data
%managagement system

%To dos:
%sort out the github/kraken situation

%% MATLAB ANALYSIS CODE LOAD IN BELOW
analysiscode_path = fullfile('/Users','andytekriwal','Desktop','git','emuanalyses_personalrepo', 'pdil');
addpath(genpath(analysiscode_path))

%% basic startup stuffs
loadOldBrkpoints = 0; %if actively working on, then set this to 0 but if opening fresh set to 1
% % %AT 1/8/20; to save out the breakpoints, use below code:
% %don't forget to manually change name of buggybrkpnt file
% b = dbstatus('-completenames');
% mainBox_AT = fullfile('/Users','andytekriwal','Desktop','git','emuanalyses_personalrepo', 'pdil','Matlab_breakpts');
% cd(mainBox_AT);
% save buggybrkpnts_PDilSandboxV2 b;
% to load in, use below:

if loadOldBrkpoints == 1
    load buggybrkpnts_PDilSandboxV2 b
    dbstop(b)
end

ptPlaceholder = 'pt_01';
behaviorFolderName = 'pt_01_postOpDay2_PDil';
behaviorFileName = 'pt_01_blockNum_1_computerTT_coop_taskoutput';
ephysfilename = 'PO_Day_04.nwb';

%% BEHAVIOR DATA LOAD IN BELOW

maindir = fullfile('/Users','andytekriwal','Box','EMU','STUDY_PDIL',ptPlaceholder, 'Behavioral Data', behaviorFolderName); % maindir = uigetdir();
cd(maindir);
load(behaviorFileName); %all behavioral structs should be loaded in as 'taskoutput'

%% Behavioral analysis


%Task navigation
for aa = 5
    timeTillStartQueries = taskoutput.timing.Time_scrn_flip_trials1(aa, 2)+taskoutput.timing.Time_scrn_flip_trials2(aa, 2)+taskoutput.timing.Time_scrn_flip_trials3(aa, 2)+taskoutput.timing.Time_scrn_flip_trials4(aa, 2)+taskoutput.timing.Time_scrn_flip_trials5(aa, 2)+taskoutput.timing.Time_postscrn_flip_trials1(aa,1)+taskoutput.timing.Time_postscrn_flip_trials2(aa,1)+taskoutput.timing.Time_postscrn_flip_trials3(aa,1)+taskoutput.timing.Time_postscrn_flip_trials4(aa,1)+taskoutput.timing.Time_postscrn_flip_trials5(aa,1);
end



timeTillEndQuery1 = taskoutput.timing.Time_query1(1) + taskoutput.timing.Time_query1_rxntime(1) + taskoutput.timing.Time_query1_Postrxntime(1);


timeTillEndQuery2 = taskoutput.timing.Time_query2(1) + taskoutput.timing.Time_query2_rxntime(1) + taskoutput.timing.Time_query2_Postrxntime(1);


sumforexample = timeTillStartQueries+timeTillEndQuery1+timeTillEndQuery2;


%6/21/20 so the above is giving us the start and stop times for when
%query1 and query2 are posed, following the 5th trial.







































%%
%%
%%
%% EPHYS DATA LOAD IN BELOW

%% add some paths
gitgenpersonalToolbox_path = fullfile('/Users','andytekriwal','Desktop','git','personalToolbox');
addpath(genpath(gitgenpersonalToolbox_path))

nwbpackage_analysis_path = fullfile('/Users','andytekriwal','Desktop','git','matnwb');
addpath(genpath(nwbpackage_analysis_path))

% Locate nwb file
maindir = fullfile('/Users','andytekriwal','Box','EMU','STUDY_PDIL',ptPlaceholder, 'SEEG', 'processed'); % maindir = uigetdir();
cd(maindir);
% load in nwb file
generateCore()
nwb_pdil = nwbRead(fullfile('/Users','andytekriwal','Box','EMU','STUDY_PDIL',ptPlaceholder, 'SEEG', 'processed', ephysfilename));
%below is helpful for navigating through the nwb files
util.nwbTree(nwb);
%% note EC has put the data under 'processing',

dat = nwb.acquisition.get('wire_14_electrode_5').data.load;
trials = nwb.intervals_trials.loadAll;

size(dat)
figure()
plot(dat(1:10000))






















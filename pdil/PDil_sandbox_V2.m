%AT 11/10/21, beginning to look at PDil analysis 

%To dos:
%sort out the github/kraken situation


%% basic startup stuffs
loadOldBrkpoints = 0; %if actively working on, then set this to 0 but if opening fresh set to 1
% % %AT 1/8/20; to save out the breakpoints, use below code:
% %don't forget to manually change name of buggybrkpnt file
% b = dbstatus('-completenames');
% mainBox_AT = fullfile('/Users','andytekriwal','Box','EMU','Matlab_breakpts');
% cd(mainBox_AT);
% save buggybrkpnts_PDilSandboxV2 b;
% to load in, use below:
if loadOldBrkpoints == 1
    load buggybrkpnts_PDilSandboxV2 b
    dbstop(b)
end


ptPlaceholder = 'pt_01';
behaviorfilename = 'pt_01_postOpDay2_PDil';
ephysfilename = 'PO_Day_04.nwb';

%% MATLAB ANALYSIS CODE LOAD IN BELOW
fullfile('/Users','andytekriwal','Desktop','git','emuanalyses_personalrepo', 'pdil');
addpath(genpath(gitgenpersonalToolbox_path))

%% BEHAVIOR DATA LOAD IN BELOW

maindir = fullfile('/Users','andytekriwal','Box','EMU','STUDY_PDIL',ptPlaceholder, 'Behavioral Data', 'processed'); % maindir = uigetdir();
cd(maindir);

%% EPHYS DATA LOAD IN BELOW
%% add some paths
gitgenpersonalToolbox_path = fullfile('/Users','andytekriwal','Desktop','git','personalToolbox');
addpath(genpath(gitgenpersonalToolbox_path))

EMUanalysis_path = fullfile('/Users','andytekriwal','Desktop','emu_analysis');
addpath(genpath(EMUanalysis_path))

% Locate nwb file 
maindir = fullfile('/Users','andytekriwal','Box','EMU','STUDY_PDIL',ptPlaceholder, 'SEEG', 'processed'); % maindir = uigetdir();
cd(maindir);
% load in nwb file
generateCore()
nwb_pdil = nwbRead(fullfile('/Users','andytekriwal','Box','EMU','STUDY_PDIL',ptPlaceholder, 'SEEG', 'processed', ephysfilename));

%% note EC has put the data under 'processing', 
dat = nwb.acquisition.get('wire_14_electrode_5').data.load;
trials = nwb.intervals_trials.loadAll;

size(dat)
figure()
plot(dat(1:10000))













    

    

    




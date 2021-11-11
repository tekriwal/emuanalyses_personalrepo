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

%% add some paths
gitgenpersonalToolbox_path = fullfile('/Users','andytekriwal','Desktop','git','personalToolbox');
addpath(genpath(gitgenpersonalToolbox_path))

EMUanalysis_path = fullfile('/Users','andytekriwal','Desktop','emu_analysis');
addpath(genpath(EMUanalysis_path))

%% load in nwb file

% Locate nwb file 
filename = 'PO_Day_04.nwb';
ptPlaceholder = 'pt_01';
maindir = fullfile('/Users','andytekriwal','Box','EMU','STUDY_PDIL',ptPlaceholder, 'SEEG', 'processed'); % maindir = uigetdir();
cd(maindir);

generateCore()

nwb_pdil = nwbRead(fullfile('/Users','andytekriwal','Box','EMU','STUDY_PDIL',ptPlaceholder, 'SEEG', 'processed', filename));

%note EC has put the data under 'processing', 
dat = nwb.acquisition.get('wire_14_electrode_5').data.load;
trials = nwb.intervals_trials.loadAll;

size(dat)
figure()
plot(dat(1:10000))













    

    




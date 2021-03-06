% https://neurodatawithoutborders.github.io/matnwb/tutorials/html/ecephys.html
% ^good tutorial

%AT 11/10/21, beginning to look at med analysis 

%To dos:
%sort out the github/kraken situation


%% basic startup stuffs
loadOldBrkpoints = 0; %if actively working on, then set this to 0 but if opening fresh set to 1
% % %AT 1/8/20; to save out the breakpoints, use below code:
% %don't forget to manually change name of buggybrkpnt file
% b = dbstatus('-completenames');
% mainBox_AT = fullfile('/Users','andytekriwal','Box','EMU','Matlab_breakpts');
% cd(mainBox_AT);
% save buggybrkpnts_medstudySandboxV2 b;
% to load in, use below:
if loadOldBrkpoints == 1
    load buggybrkpnts_medstudySandboxV2 b
    dbstop(b)
end


%% add some paths
gitgenpersonalToolbox_path = fullfile('/Users','andytekriwal','Desktop','git','personalToolbox');
addpath(genpath(gitgenpersonalToolbox_path))

EMUanalysis_path = fullfile('/Users','andytekriwal','Desktop','emu_analysis');
addpath(genpath(EMUanalysis_path))

%% load in nwb file

% Locate nwb file 
filename = 'patient2535.nwb';
ptPlaceholder = 'pt_2535';
maindir = fullfile('/Users','andytekriwal','Box','EMU','STUDY_MEDS', ptPlaceholder); % maindir = uigetdir();
cd(maindir);

generateCore()

nwb = nwbRead( fullfile('/Users','andytekriwal','Box','EMU','STUDY_MEDS', ptPlaceholder, filename) ); % maindir = uigetdir();

util.nwbTree(nwb);

nwb.processing

nwb.processing.get('ecephys')

nwb.processing.get('ecephys').nwbdatainterface.get('Lpar')

nwb.processing.get('ecephys').nwbdatainterface.get('Lpar').data.load

%














    
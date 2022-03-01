

%starting on 2/3/22

%README: purpose of script is to generate behavioral comparisons between
%the six task participants we have as of 2/3/22. There are a few key places
%where we can alter the code. Namely, which trials in the block we want to
%investigate ('endtrial' and 'starttrial') and if we want to
%skip specific blocks from a patient (specifically pt 4, but atm we've
%elected to keep all data from pt)




dbstop if error
%below two vars tell which trials in a given session we should concentrate
%on, may want to change based on exactly what we're looking at
starttrial = 1;
endtrial = 15;


pdilBehaviorStruct = struct();

%Below lists what are possible vars to input given how things are labeled
%pt_01; day = 1 or 2, with blocks 1-2 on day 1, and blocks 3-6 on day 2
%pt_02; day = 1 or 2, with blocks 1-2 on day 1, and blocks 3-4 on day 2
%pt_03; day = 1, with blocks 1-4

%labeled as 'pt#' 'day#' and then 'block#'
%pt_01
inputlist(1,:) = [1;1;1];
inputlist(2,:) = [1;1;2];
inputlist(3,:) = [1;2;3];
inputlist(4,:) = [1;2;4];
inputlist(5,:) = [1;2;5];
inputlist(6,:) = [1;2;6];
%pt_02
inputlist(7,:) = [2;1;1];
inputlist(8,:) = [2;1;2];
inputlist(9,:) = [2;2;3];
inputlist(10,:) = [2;2;4];
%pt_03
inputlist(11,:) = [3;1;1];
inputlist(12,:) = [3;1;2];
inputlist(13,:) = [3;1;3];
inputlist(14,:) = [3;1;4];
%pt_04
inputlist(15,:) = [4;1;1];
inputlist(16,:) = [4;1;2];
inputlist(17,:) = [4;1;3];
%pt_05
inputlist(18,:) = [5;1;1];
inputlist(19,:) = [5;2;2];
%pt_06
inputlist(20,:) = [6;1;1];%note that I'm not sure if these blocks were indeed on the same day, but pretty sure
inputlist(21,:) = [6;1;2];
inputlist(22,:) = [6;1;3];
inputlist(23,:) = [6;1;4];


for blocknumberindex = 1:length(inputlist)
    
    keyinputlist = inputlist(blocknumberindex,:);
    
    ptPlaceholder = strcat('pt_0', mat2str(keyinputlist(1)));
    day = keyinputlist(2);
    block = keyinputlist(3);
    
    identifier = strcat(ptPlaceholder, '__day', mat2str(day), '__block', mat2str(block));
    
    pdilBehaviorStruct(blocknumberindex).id = identifier;
    %%
    
    %add key repos
    analysiscode_path = fullfile('/Users','andytekriwal','Desktop','git','emuanalyses_personalrepo', 'pdil');
    addpath(genpath(analysiscode_path))
    
    gitgenpersonalToolbox_path = fullfile('/Users','andytekriwal','Desktop','git','personalToolbox');
    addpath(genpath(gitgenpersonalToolbox_path))
    
    
    % to save out the breakpoints, use below code
    saveBrkpoints = 0; %if actively working on, then set this to 0 but if opening fresh set to 1
    if saveBrkpoints == 1
        b = dbstatus('-completenames');
        mainBox_AT = fullfile('/Users','andytekriwal','Desktop','git','emuanalyses_personalrepo', 'pdil','Matlab_breakpts');
        cd(mainBox_AT);
        delete buggybrkpnts_pdilBehaviorAnalysis_V1.mat
        save buggybrkpnts_pdilBehaviorAnalysis_V1.mat b;
    end
    
    loadOldBrkpoints = 0; %if actively working on, then set this to 0 but if opening fresh set to 1
    if loadOldBrkpoints == 1
        load buggybrkpnts_pdilBehaviorAnalysis_V1.mat b
        dbstop(b)
    end
    
    
    
    %%
    %identify behavior file names and
    if strcmp(ptPlaceholder,'pt_01')
        if day == 1
            behaviorFolderName = 'pt_01_postOpDay2_PDil';
            if block == 1
                behaviorFileName = 'pt_01_blockNum_1_computerTT_coop_taskoutput';
            elseif block == 2
                behaviorFileName = 'pt_01_blockNum_2_humanTT_coop_taskoutput';
            end
        elseif day == 2
            behaviorFolderName = 'pt_01_postOpDay4_PDil';
            if block == 3
                behaviorFileName = 'pt_01_blockNum_3_computerTT_defect_taskoutput';
            elseif block == 4
                behaviorFileName = 'pt_01_blockNum_4_humanTT_defect_taskoutput';
            elseif block == 5
                behaviorFileName = 'pt_01_blockNum_5_humanTT_coop_taskoutput';
            elseif block == 6
                behaviorFileName = 'pt_01_blockNum_6_computerTT_coop_taskoutput';
            end
        end
    elseif strcmp(ptPlaceholder,'pt_02')
        if day == 1
            behaviorFolderName = 'pt_02_postOpDay4_PDil';
            if block == 1
                behaviorFileName = 'pt_02_blockNum_1_computerTT_coop_taskoutput';
            elseif block == 2
                behaviorFileName = 'pt_02_blockNum_2_humanTT_coop_taskoutput';
            end
        elseif day == 2
            behaviorFolderName = 'pt_02_postOpDay5_PDil';
            if block == 3
                behaviorFileName = 'pt_02_blockNum_3_computerTT_defect_taskoutput';
            elseif block == 4
                behaviorFileName = 'pt_02_blockNum_4_humanTT_defect_taskoutput';
            end
        end
    elseif strcmp(ptPlaceholder,'pt_03')
        if day == 1
            behaviorFolderName = 'pt_03_postOpDay2_PDil';
            if block == 1
                behaviorFileName = 'pt_03_blockNum_1_computerTT_coop_taskoutput';
            elseif block == 2
                behaviorFileName = 'pt_03_blockNum_2_humanTT_coop_taskoutput';
            elseif block == 3
                behaviorFileName = 'pt_03_blockNum_3_computerTT_defect_taskoutput';
            elseif block == 4
                behaviorFileName = 'pt_03_blockNum_4_humanTT_defect_taskoutput';
            end
        end
    elseif strcmp(ptPlaceholder,'pt_04')
        if day == 1
            behaviorFolderName = 'pt_04_postOpDay2_PDil';
            if block == 1
                behaviorFileName = 'pt_04_blockNum_1_computerTT_coop_taskoutput';
            elseif block == 2
                behaviorFileName = 'pt_04_blockNum_2_humanTT_coop_taskoutput';
            elseif block == 3
                behaviorFileName = 'pt_04_blockNum_3_computerTT_defect_taskoutput';
            end
        end
    elseif strcmp(ptPlaceholder,'pt_05')
        if day == 1
            behaviorFolderName = 'pt_05_postOpDay2_PDil';
            if block == 1
                behaviorFileName = 'pt_05_blockNum_1_humanTT_coop_taskoutput';
            end
        elseif day == 2
            behaviorFolderName = 'pt_05_postOpDay3_PDil';
            if block == 2
                behaviorFileName = 'pt_05_blockNum_2_computerTT_coop_taskoutput';
            end
        end
    elseif strcmp(ptPlaceholder,'pt_06')
        if day == 1
            behaviorFolderName = 'pt_06_postOpDay3_PDil';
            if block == 1
                behaviorFileName = 'pt_06_blockNum_1_humanTT_coop_taskoutput';
            elseif block == 2
                behaviorFileName = 'pt_06_blockNum_2_computerTT_coop_taskoutput';
            elseif block == 3
                behaviorFileName = 'pt_06_blockNum_3_humanTT_defect_taskoutput';
            elseif block == 4
                behaviorFileName = 'pt_06_blockNum_4_computerTT_defect_taskoutput';
            end
        end
    end
    
    
    pdilBehaviorStruct(blocknumberindex).id_behaviorFileName = behaviorFileName;
    
    if strcmp(behaviorFileName(18:24), 'humanTT')
        pdilBehaviorStruct(blocknumberindex).opponent = "human";
    elseif strcmp(behaviorFileName(18:24), 'compute')
        pdilBehaviorStruct(blocknumberindex).opponent = "computer";
    end
    
    if strcmp(behaviorFileName(27:32), 'T_coop') || strcmp(behaviorFileName(27:32), 'oop_ta')
        pdilBehaviorStruct(blocknumberindex).opponentstrat = "coop";
    elseif strcmp(behaviorFileName(27:32), 'T_defe') || strcmp(behaviorFileName(27:32), 'efect_')
        pdilBehaviorStruct(blocknumberindex).opponentstrat = "defect";
    end
    
    %load behavioral files
    maindir = fullfile('/Users','andytekriwal','Library','CloudStorage','Box-Box', 'EMU','STUDY_PDIL',ptPlaceholder, 'Behavioral Data', behaviorFolderName); % maindir = uigetdir();
    cd(maindir);
    load(behaviorFileName); %all behavioral structs should be loaded in as 'taskoutput'
    
    % if taskoutput.choice_matrix_player1(or2)_logical == 1, means cooperate; whereas == 2 means defect
    pdilBehaviorStruct(blocknumberindex).choice_matrix_player1_logical = taskoutput.choice_matrix_player1_logical;
    pdilBehaviorStruct(blocknumberindex).choice_matrix_player2_logical = taskoutput.choice_matrix_player2_logical;
    
    %Below are output vars for the coop/trust questions that occur after every 5 blocks
    pdilBehaviorStruct(blocknumberindex).probeResponse_playerAcoop = taskoutput.probeResponse_playerAcoop;
    pdilBehaviorStruct(blocknumberindex).probeResponse_playerAtrust = taskoutput.probeResponse_playerAtrust;
    
    
    
    
    
    blocknumber = str2double(behaviorFileName(16));
    %sanity check regarding block# as extracted from behavioral file and block#
    %as indicated in key vars at top of script
    if blocknumber ~= block
        disp('ERROR: mismatch in stated block numb and extracted block numb')
        2+2
    end
    
    %%
    %QUESTIONS EPOCH
    taskoutput.timing_sumTictocs_trial_n; %keep in mind this is the output var for the sum of tictocs in task output
    
    
    for aa = 1:15
        typicalTrialEventsDuration(aa) = taskoutput.timing.Time_scrn_flip_trials1(aa, 2) + taskoutput.timing.Time_scrn_flip_trials2(aa, 2) + taskoutput.timing.Time_scrn_flip_trials3(aa, 2) + taskoutput.timing.Time_scrn_flip_trials4(aa, 2) + taskoutput.timing.Time_scrn_flip_trials5(aa, 2) + taskoutput.timing.Time_postscrn_flip_trials1(aa,1) + taskoutput.timing.Time_postscrn_flip_trials2(aa,1) + taskoutput.timing.Time_postscrn_flip_trials3(aa,1) + taskoutput.timing.Time_postscrn_flip_trials4(aa,1) + taskoutput.timing.Time_postscrn_flip_trials5(aa,1);
    end
    
    pdilBehaviorStruct(blocknumberindex).typicalTrialEventsDuration = typicalTrialEventsDuration;
    
    %rxn time for query1 (how cooperative have you been) across the block, as well as query2 (how trustworth is
    %opp)
    
    pdilBehaviorStruct(blocknumberindex).Time_query1_rxntime = taskoutput.timing.Time_query1_rxntime;
    pdilBehaviorStruct(blocknumberindex).Time_query2_rxntime = taskoutput.timing.Time_query2_rxntime;
    
    initiateTrial = zeros(1,15);CandDpromptAppears = zeros(1,15);CandDpromptResponseRequest = zeros(1,15);CandDpromptResponseGiven = zeros(1,15);spacebarToSeeOpponentsChoice = zeros(1,15);opponentsChoiceOnScreen = zeros(1,15);doneSeeingOpponentsChoice = zeros(1,15);
    
    for aa = 1:15
        initiateTrial(aa) = taskoutput.timing.Time_scrn_flip_trials1(aa,2)+taskoutput.timing.Time_postscrn_flip_trials1(aa);
        CandDpromptAppears(aa) = initiateTrial(aa) + taskoutput.timing.Time_scrn_flip_trials2(aa,1);
        CandDpromptResponseRequest(aa) = initiateTrial(aa) + taskoutput.timing.Time_scrn_flip_trials2(aa,2); %this should be a second later than the above
        CandDpromptResponseGiven(aa) = CandDpromptResponseRequest(aa) + taskoutput.timing.Time_postscrn_flip_trials2(aa);
        spacebarToSeeOpponentsChoice(aa) = CandDpromptResponseGiven(aa)+ taskoutput.timing.Time_scrn_flip_trials3(aa,2)+taskoutput.timing.Time_postscrn_flip_trials3(aa)+taskoutput.timing.Time_scrn_flip_trials4(aa,2);
        opponentsChoiceOnScreen(aa) = spacebarToSeeOpponentsChoice(aa) + taskoutput.timing.Time_postscrn_flip_trials4(aa)+taskoutput.timing.Time_scrn_flip_trials5(aa,1);
        doneSeeingOpponentsChoice(aa) = spacebarToSeeOpponentsChoice(aa) + taskoutput.timing.Time_postscrn_flip_trials4(aa)+taskoutput.timing.Time_scrn_flip_trials5(aa,2);
    end
    
    pdilBehaviorStruct(blocknumberindex).trialStart_to_CorDselection_startstop = ... %time from when trial was initiated with spacebar and when the either C or D response was given by the study pariticpant
        initiateTrial-CandDpromptResponseGiven;
    
    pdilBehaviorStruct(blocknumberindex).CorDselection_to_doneSeeingOpponentsChoice_startstop = ... %time from when C/D response given and until clicked past opponent's response
        doneSeeingOpponentsChoice - CandDpromptResponseGiven;
    
    pdilBehaviorStruct(blocknumberindex).CorDrxntime = ... %reaction time for participant
        CandDpromptResponseGiven - CandDpromptResponseRequest;
    
end

Player1_percCooperated = zeros(length(pdilBehaviorStruct),1);
Player2_percCooperated = Player1_percCooperated; 
CorDrxntime = Player1_percCooperated; 
CorDselection_to_doneSeeingOpponentsChoice_startstop = Player1_percCooperated; 

totaltrials = endtrial - starttrial + 1;
for mm = 1:length(pdilBehaviorStruct)
    
    Player1_percCooperated(mm) = ((sum(pdilBehaviorStruct(mm).choice_matrix_player1_logical(starttrial:endtrial)))/totaltrials)*100;
    Player2_percCooperated(mm) = ((sum(pdilBehaviorStruct(mm).choice_matrix_player2_logical(starttrial:endtrial)))/totaltrials)*100;
%     CorDrxntime(mm) = ((sum(pdilBehaviorStruct(mm).CorDrxntime(starttrial:endtrial)))/totaltrials);
    CorDrxntime(mm) = median(pdilBehaviorStruct(mm).CorDrxntime(starttrial:endtrial));
%     CorDselection_to_doneSeeingOpponentsChoice_startstop(mm) = ((sum(pdilBehaviorStruct(mm).CorDselection_to_doneSeeingOpponentsChoice_startstop(starttrial:endtrial)))/totaltrials);
    CorDselection_to_doneSeeingOpponentsChoice_startstop(mm) = median(pdilBehaviorStruct(mm).CorDselection_to_doneSeeingOpponentsChoice_startstop(starttrial:endtrial));
  
    
    opponentstrat(mm,1) = pdilBehaviorStruct(mm).opponentstrat;
    opponent(mm,1) = pdilBehaviorStruct(mm).opponent;
    
    probeResponse_playerAcoop(1,mm) = str2num(pdilBehaviorStruct(mm).probeResponse_playerAcoop(1));
    probeResponse_playerAcoop(2,mm) = str2num(pdilBehaviorStruct(mm).probeResponse_playerAcoop(2));
    probeResponse_playerAcoop(3,mm) = str2num(pdilBehaviorStruct(mm).probeResponse_playerAcoop(3));
    probeResponse_playerAcoop(4, mm) = (str2num(pdilBehaviorStruct(mm).probeResponse_playerAcoop(1))+str2num(pdilBehaviorStruct(mm).probeResponse_playerAcoop(2))+str2num(pdilBehaviorStruct(mm).probeResponse_playerAcoop(3)))/3;

     probeResponse_playerAtrust(1,mm) = str2num(pdilBehaviorStruct(mm).probeResponse_playerAtrust(1));
    probeResponse_playerAtrust(2,mm) = str2num(pdilBehaviorStruct(mm).probeResponse_playerAtrust(2));
    probeResponse_playerAtrust(3,mm) = str2num(pdilBehaviorStruct(mm).probeResponse_playerAtrust(3));
    probeResponse_playerAtrust(4, mm) = (str2num(pdilBehaviorStruct(mm).probeResponse_playerAtrust(1))+str2num(pdilBehaviorStruct(mm).probeResponse_playerAtrust(2))+str2num(pdilBehaviorStruct(mm).probeResponse_playerAtrust(3)))/3;
        
    Slope_probeResponse_playerAcoop(1,mm) = str2num(pdilBehaviorStruct(mm).probeResponse_playerAcoop(2)) - str2num(pdilBehaviorStruct(mm).probeResponse_playerAcoop(1));       
    Slope_probeResponse_playerAcoop(2,mm) = str2num(pdilBehaviorStruct(mm).probeResponse_playerAcoop(3)) - str2num(pdilBehaviorStruct(mm).probeResponse_playerAcoop(1));       
    Slope_probeResponse_playerAcoop(3,mm) = str2num(pdilBehaviorStruct(mm).probeResponse_playerAcoop(3)) - str2num(pdilBehaviorStruct(mm).probeResponse_playerAcoop(2));       
        
    Slope_probeResponse_playerAtrust(1,mm) = str2num(pdilBehaviorStruct(mm).probeResponse_playerAtrust(2)) - str2num(pdilBehaviorStruct(mm).probeResponse_playerAtrust(1));       
    Slope_probeResponse_playerAtrust(2,mm) = str2num(pdilBehaviorStruct(mm).probeResponse_playerAtrust(3)) - str2num(pdilBehaviorStruct(mm).probeResponse_playerAtrust(1));       
    Slope_probeResponse_playerAtrust(3,mm) = str2num(pdilBehaviorStruct(mm).probeResponse_playerAtrust(3)) - str2num(pdilBehaviorStruct(mm).probeResponse_playerAtrust(2));       

% AT below 2/14/22. more sophisticated method of determining coop/defects
    %fraction of opportunities that 'revenge' was taken
    matrix1 = pdilBehaviorStruct(mm).choice_matrix_player1_logical(starttrial:endtrial)
    matrix2 = pdilBehaviorStruct(mm).choice_matrix_player2_logical(starttrial:endtrial)
    %find where opponent defected (but don't look at last trial)
    finddef = find(matrix2(1:(endtrial-1))==0)
    subjectresponseIndex = finddef + 1;
    findresponses_binary = matrix1(subjectresponseIndex);
    timesdefected_aftedDefectedOn = length(find(findresponses_binary==0))
    fraction_timesdefected_aftedDefectedOn(mm) = timesdefected_aftedDefectedOn/length(finddef)
 
    %fraction of times defection was unprovoked
    matrix1 = pdilBehaviorStruct(mm).choice_matrix_player1_logical(starttrial:endtrial)
    matrix2 = pdilBehaviorStruct(mm).choice_matrix_player2_logical(starttrial:endtrial)
    %find where subject defected (not counting first trial)
    finddefs = find(matrix1(2:end)==0) %so in reality, this pulls out the trialindex for the trial BEFORE when subject defects
    oppresponseIndex = finddefs;
    findresponses_binary2 = matrix2(oppresponseIndex);
    timesdefected_withoutcause = length(find(findresponses_binary2==1))
    fraction_timesdefected_withoutcause(mm) = timesdefected_withoutcause/length(finddefs)
    
    
    if length(finddefs) == 0
        fraction_timesdefected_withoutcause(mm) = 0;
    end
    
    
    %find where opponent defected
    
%     for zz = 1:(length(matrix1)-1)
%         if matrix2(zz) == 0
%         end
%     end
      
end

%%
szz = 600;
longersideInput = 600; %9/4/20, AT's attempt at Golden ratio calcs for size


skippt05 = 0; %set to 0 to skip no sessions, set to .5 to skip the first session of pt 5 (the one he got competency query wrong), set to 1 to skip all of pt 5.
% AT, as of 2/14/22, we want to NOT skip any sessions, so keep set to 0

if skippt05 == 1
    skipPT05 = 3;
    removept4Index = [17, 16, 15];
    for i = 1:length(removept4Index)
        
        opponentstrat(removept4Index(i)) = [];
        opponent(removept4Index(i)) = [];
        Player1_percCooperated(removept4Index(i)) = [];
        Player2_percCooperated(removept4Index(i)) = [];
        CorDrxntime(removept4Index(i)) = [];
        Slope_probeResponse_playerAcoop(:, removept4Index(i)) = [];
        Slope_probeResponse_playerAtrust(:, removept4Index(i)) = [];
    end
    
elseif skippt05 == .5
     skipPT05 = 1;
    removept4Index = 15;
    for i = 1:length(removept4Index)
        
        opponentstrat(removept4Index(i)) = [];
        opponent(removept4Index(i)) = [];
        Player1_percCooperated(removept4Index(i)) = [];
        Player2_percCooperated(removept4Index(i)) = [];
        CorDrxntime(removept4Index(i)) = [];
        Slope_probeResponse_playerAcoop(:, removept4Index(i)) = [];
        Slope_probeResponse_playerAtrust(:, removept4Index(i)) = [];

    end
       
elseif skippt05 == 0
    skipPT05 = 0;
%     removept4Index = 15;
%     for i = 1:length(removept4Index)
%         
%         opponentstrat(removept4Index(i)) = [];
%         opponent(removept4Index(i)) = [];
%         Player1_percCooperated(removept4Index(i)) = [];
%         Player2_percCooperated(removept4Index(i)) = [];
%         CorDrxntime(removept4Index(i)) = [];
%     end
%     
    
end

%%
%%
%%
%% PLOTTING!

%% BELOW IS LOOKING AT THE PERCENT OF TRIALS THE SUBJECT COOPERATED ON
Ylim = [0 100];
keymetric = Player1_percCooperated;
%keymetric = probeResponse_playerAcoop 
%keymetric = probeResponse_playerAtrust
ylabeltext = '% of trials in a block the player cooperated'; %, 'FontSize', 14);
Ylim = [0 100];

keymetric = fraction_timesdefected_aftedDefectedOn*100
ylabeltext = '% of opportunities revenge taken'; %, 'FontSize', 14);

keymetric = fraction_timesdefected_withoutcause*100 %this is very interested for vsCoop and vsDefect
ylabeltext = '% of defects that were unprovoked'; %, 'FontSize', 14);


% % Below compares the play of subject and opponents (not too useful)
% keymetric = keymetric;
% input2 = Player2_percCooperated;
% xlabel1text = 'Subject'; 
% xlabel2text = 'Opponent';
% pairedplotting_helperfx(keymetric, input2, xlabel1text, xlabel2text, ylabeltext)


% below is plot comparing cooperates to defects for whatever key variable
index = find(opponentstrat=="coop");
input1 = keymetric(index);
index = find(opponentstrat=="defect");
input2 = keymetric(index);
xlabel1text = 'vsCoop'; 
xlabel2text = 'vsDefect';
ylabeltext = ylabeltext; %, 'FontSize', 14);
notpairedplotting_helperfx(input1, input2, xlabel1text, xlabel2text, ylabeltext, Ylim)

%below is for comparing within subject for against coop or defect
%strategies

%within subject vscoop VS vsdefect
keymetric = keymetric; clear pairedmat
xlabel1text = 'vsCoop'; xlabel2text = 'vsDefect';
ylabeltext = ylabeltext; %, 'FontSize', 14);
%pt1 below
pairedmat(1,1) = mean([keymetric(1), keymetric(2), keymetric(5), keymetric(6)]);%vscoop 
pairedmat(1,2) = mean([keymetric(3), keymetric(4)]);%vsdefect
%pt2 below
pairedmat(2,1) = mean([keymetric(7), keymetric(8)]);%vscoop
pairedmat(2,2) = mean([keymetric(9), keymetric(10)]);%vsdefect
%pt3 below
pairedmat(3,1) = mean([keymetric(11), keymetric(12)]);%vscoop 
pairedmat(3,2) = mean([keymetric(13), keymetric(14)]);%vsdefect
% %pt5 below %has no blocks against a defective strat
% pairedmat(4,1) = mean([keymetric(18), keymetric(19)]);%vscoop (no pt 4 here)
% pairedmat(4,2) = nan;%vsdefect - none such trials 
%pt6 below
pairedmat(4,1) = mean([keymetric(20-skipPT05), keymetric(21-skipPT05)]);%vscoop 
pairedmat(4,2) = mean([keymetric(22-skipPT05), keymetric(23-skipPT05)]);%vsdefect 
%AT adding below 2/11/22 to include data from pt4
if skippt05 == 0
    pairedmat(5,1) = mean([keymetric(15-skipPT05), keymetric(16-skipPT05)]);%vscoop 
    pairedmat(5,2) = keymetric(17-skipPT05);%vsdefect 
elseif skippt05 == 0.5
    pairedmat(5,1) = keymetric(16-skipPT05);%vscoop 
    pairedmat(5,2) = keymetric(17-skipPT05);%vsdefect 
end

input1 = pairedmat(:,1); input2 = pairedmat(:,2);
pairedplotting_helperfx(input1, input2, xlabel1text, xlabel2text, ylabeltext, Ylim)







% below is plot comparing cooperates to defects for whatever key variable
index = find(opponent=="human");
input1 = keymetric(index);
index = find(opponent=="computer");
input2 = keymetric(index);
xlabel1text = 'vsHuman'; 
xlabel2text = 'vsComp';
ylabeltext = ylabeltext; %, 'FontSize', 14);
notpairedplotting_helperfx(input1, input2, xlabel1text, xlabel2text, ylabeltext, Ylim)

%within subject vsperson VS vscomputer
keymetric = keymetric; clear pairedmat
xlabel1text = 'vsHuman'; xlabel2text = 'vsComp';
ylabeltext = ylabeltext; %, 'FontSize', 14);
%pt1 below
pairedmat(1,1) = mean([keymetric(2), keymetric(4), keymetric(5)]);%vsperson 
pairedmat(1,2) = mean([keymetric(1), keymetric(3), keymetric(6)]);%vscomp
%pt2 below
pairedmat(2,1) = mean([keymetric(8), keymetric(10)]);%vsperson
pairedmat(2,2) = mean([keymetric(7), keymetric(9)]);%vscomp
%pt3 below
pairedmat(3,1) = mean([keymetric(12), keymetric(14)]);%vsperson 
pairedmat(3,2) = mean([keymetric(11), keymetric(13)]);%vscomp
%pt5 below
pairedmat(4,1) = keymetric(18-skipPT05);%vsperson (remember, no pt 4 here)
pairedmat(4,2) = keymetric(19-skipPT05);%vscomputer (remember, no pt 4 here)
%pt6 below
pairedmat(5,1) = mean([keymetric(20-skipPT05), keymetric(22-skipPT05)]);%vsperson (remember, no pt 4 here)
pairedmat(5,2) = mean([keymetric(21-skipPT05), keymetric(23-skipPT05)]);%vscomputer (remember, no pt 4 here)
%AT adding below 2/11/22 to include data from pt4
if skippt05 == 0
    pairedmat(6,1) = keymetric(16-skipPT05);%vshuman 
    pairedmat(6,2) = mean([keymetric(15-skipPT05), keymetric(17-skipPT05)]);%vscomputer 
elseif skippt05 == 0.5
    pairedmat(6,1) = keymetric(16-skipPT05);%vshuman 
    pairedmat(6,2) = keymetric(17-skipPT05);%vscomputer 
end



input1 = pairedmat(:,1); input2 = pairedmat(:,2);
pairedplotting_helperfx(input1, input2, xlabel1text, xlabel2text, ylabeltext, Ylim)





2+2


%% BELOW IS LOOKING AT RXN TIMES

keymetric = CorDrxntime;
ylabeltext = 'Time for C/D selection (seconds)'; %, 'FontSize', 14);
Ylim = [0 6];

keymetric = CorDselection_to_doneSeeingOpponentsChoice_startstop;
ylabeltext = 'Time spent awaiting opp selection (seconds)'; %, 'FontSize', 14);
Ylim = [6 10];



% % Below compares the play of subject and opponents (not too useful)
% keymetric = keymetric;
% input2 = Player2_percCooperated;
% xlabel1text = 'Subject'; 
% xlabel2text = 'Opponent';
% pairedplotting_helperfx(keymetric, input2, xlabel1text, xlabel2text, ylabeltext)


% below is plot comparing cooperates to defects for whatever key variable
index = find(opponentstrat=="coop");
input1 = keymetric(index);
index = find(opponentstrat=="defect");
input2 = keymetric(index);
xlabel1text = 'vsCoop'; 
xlabel2text = 'vsDefect';
ylabeltext = ylabeltext; %, 'FontSize', 14);
notpairedplotting_helperfx(input1, input2, xlabel1text, xlabel2text, ylabeltext, Ylim)

%below is for comparing within subject for against coop or defect
%strategies

%within subject vscoop VS vsdefect
keymetric = keymetric; clear pairedmat
xlabel1text = 'vsCoop'; xlabel2text = 'vsDefect';
ylabeltext = ylabeltext; %, 'FontSize', 14);
%pt1 below
pairedmat(1,1) = mean([keymetric(1), keymetric(2), keymetric(5), keymetric(6)]);%vscoop 
pairedmat(1,2) = mean([keymetric(3), keymetric(4)]);%vsdefect
%pt2 below
pairedmat(2,1) = mean([keymetric(7), keymetric(8)]);%vscoop
pairedmat(2,2) = mean([keymetric(9), keymetric(10)]);%vsdefect
%pt3 below
pairedmat(3,1) = mean([keymetric(11), keymetric(12)]);%vscoop 
pairedmat(3,2) = mean([keymetric(13), keymetric(14)]);%vsdefect
% %pt5 below %has no blocks against a defective strat
% pairedmat(4,1) = mean([keymetric(18), keymetric(19)]);%vscoop (no pt 4 here)
% pairedmat(4,2) = nan;%vsdefect - none such trials 
%pt6 below
pairedmat(4,1) = mean([keymetric(20-skipPT05), keymetric(21-skipPT05)]);%vscoop 
pairedmat(4,2) = mean([keymetric(22-skipPT05), keymetric(23-skipPT05)]);%vsdefect 
%AT adding below 2/11/22 to include data from pt4
if skippt05 == 0
    pairedmat(5,1) = mean([keymetric(15-skipPT05), keymetric(16-skipPT05)]);%vscoop 
    pairedmat(5,2) = keymetric(17-skipPT05);%vsdefect 
elseif skippt05 == 0.5
    pairedmat(5,1) = keymetric(16-skipPT05);%vscoop 
    pairedmat(5,2) = keymetric(17-skipPT05);%vsdefect 
end

input1 = pairedmat(:,1); input2 = pairedmat(:,2);
pairedplotting_helperfx(input1, input2, xlabel1text, xlabel2text, ylabeltext, Ylim)







% below is plot comparing cooperates to defects for whatever key variable
index = find(opponent=="human");
input1 = keymetric(index);
index = find(opponent=="computer");
input2 = keymetric(index);
xlabel1text = 'vsHuman'; 
xlabel2text = 'vsComp';
ylabeltext = ylabeltext; %, 'FontSize', 14);
notpairedplotting_helperfx(input1, input2, xlabel1text, xlabel2text, ylabeltext, Ylim)

%within subject vsperson VS vscomputer
keymetric = keymetric; clear pairedmat
xlabel1text = 'vsHuman'; xlabel2text = 'vsComp';
ylabeltext = ylabeltext; %, 'FontSize', 14);
%pt1 below
pairedmat(1,1) = mean([keymetric(2), keymetric(4), keymetric(5)]);%vsperson 
pairedmat(1,2) = mean([keymetric(1), keymetric(3), keymetric(6)]);%vscomp
%pt2 below
pairedmat(2,1) = mean([keymetric(8), keymetric(10)]);%vsperson
pairedmat(2,2) = mean([keymetric(7), keymetric(9)]);%vscomp
%pt3 below
pairedmat(3,1) = mean([keymetric(12), keymetric(14)]);%vsperson 
pairedmat(3,2) = mean([keymetric(11), keymetric(13)]);%vscomp
%pt5 below
pairedmat(4,1) = keymetric(18-skipPT05);%vsperson (remember, no pt 4 here)
pairedmat(4,2) = keymetric(19-skipPT05);%vscomputer (remember, no pt 4 here)
%pt6 below
pairedmat(5,1) = mean([keymetric(20-skipPT05), keymetric(22-skipPT05)]);%vsperson (remember, no pt 4 here)
pairedmat(5,2) = mean([keymetric(21-skipPT05), keymetric(23-skipPT05)]);%vscomputer (remember, no pt 4 here)
%AT adding below 2/11/22 to include data from pt4
if skippt05 == 0
    pairedmat(6,1) = keymetric(16-skipPT05);%vshuman 
    pairedmat(6,2) = mean([keymetric(15-skipPT05), keymetric(17-skipPT05)]);%vscomputer 
elseif skippt05 == 0.5
    pairedmat(6,1) = keymetric(16-skipPT05);%vshuman 
    pairedmat(6,2) = keymetric(17-skipPT05);%vscomputer 
end



input1 = pairedmat(:,1); input2 = pairedmat(:,2);
pairedplotting_helperfx(input1, input2, xlabel1text, xlabel2text, ylabeltext, Ylim)





2+2



%% BELOW IS LOOKING AT CHANGES IN RESPONSES TO THE PROBE QUESTIONS

%below is just the probe response questions themselves:
Ylim = [1 9];
 keymetric =    probeResponse_playerAcoop(1,:); 
 ylabeltext = 'Self assessed coop level (block 1)';
 keymetric =    probeResponse_playerAcoop(2,:); 
 ylabeltext = 'Self assessed coop level (block 2)'; 
  keymetric =    probeResponse_playerAcoop(3,:); 
 ylabeltext = 'Self assessed coop level (block 3)';
  keymetric =    probeResponse_playerAcoop(4,:); 
 ylabeltext = 'Self assessed coop level (ave of all blocks)';
 
 keymetric =    probeResponse_playerAtrust(1,:); 
 ylabeltext = 'Trust in opponent (block 1)';
 keymetric =    probeResponse_playerAtrust(2,:); 
 ylabeltext = 'Trust in opponent (block 2)'; 
  keymetric =    probeResponse_playerAtrust(3,:); 
 ylabeltext = 'Trust in opponent (block 3)';
  keymetric =    probeResponse_playerAtrust(4,:); 
 ylabeltext = 'Trust in opponent (ave of all blocks)';
 
 
% below we want to calculate the 'slope' of the line that connects the baseine

Ylim = [-5 5];
% %pick from below pairs of options;
keymetric = Slope_probeResponse_playerAtrust(1,:); 
ylabeltext = 'Change in trust (block 1 to 2)'; %, 'FontSize', 14);
keymetric = Slope_probeResponse_playerAtrust(2,:); 
ylabeltext = 'Change in trust (block 1 to 3)'; %, 'FontSize', 14);
keymetric = Slope_probeResponse_playerAtrust(3,:); 
ylabeltext = 'Change in trust (block 2 to 3)'; %, 'FontSize', 14);
keymetric = mean([Slope_probeResponse_playerAtrust(1,:);Slope_probeResponse_playerAtrust(2,:)],1); 
ylabeltext = 'Ave change in trust (1 to block2, 1 to block3)'; %, 'FontSize', 14);


keymetric = Slope_probeResponse_playerAcoop(1,:); 
ylabeltext = 'Change in self assessed coop (block 1 to 2)'; %, 'FontSize', 14);
keymetric = Slope_probeResponse_playerAcoop(2,:); 
ylabeltext = 'Change in self assessed coop (block 1 to 3)'; %, 'FontSize', 14);
keymetric = Slope_probeResponse_playerAcoop(3,:); 
ylabeltext = 'Change in self assessed coop (block 2 to 3)'; %, 'FontSize', 14);


% below is plot comparing cooperates to defects for whatever key variable
index = find(opponentstrat=="coop");
input1 = keymetric(index);
index = find(opponentstrat=="defect");
input2 = keymetric(index);
xlabel1text = 'vsCoop'; 
xlabel2text = 'vsDefect';
ylabeltext = ylabeltext; %, 'FontSize', 14);
notpairedplotting_helperfx(input1, input2, xlabel1text, xlabel2text, ylabeltext, Ylim)


%below is for comparing within subject for against coop or defect
%strategies

%within subject vscoop VS vsdefect
keymetric = keymetric; clear pairedmat
xlabel1text = 'vsCoop'; xlabel2text = 'vsDefect';
ylabeltext = ylabeltext; %, 'FontSize', 14);
%pt1 below
pairedmat(1,1) = mean([keymetric(1), keymetric(2), keymetric(5), keymetric(6)]);%vscoop 
pairedmat(1,2) = mean([keymetric(3), keymetric(4)]);%vsdefect
%pt2 below
pairedmat(2,1) = mean([keymetric(7), keymetric(8)]);%vscoop
pairedmat(2,2) = mean([keymetric(9), keymetric(10)]);%vsdefect
%pt3 below
pairedmat(3,1) = mean([keymetric(11), keymetric(12)]);%vscoop 
pairedmat(3,2) = mean([keymetric(13), keymetric(14)]);%vsdefect
% %pt5 below %has no blocks against a defective strat
% pairedmat(4,1) = mean([keymetric(18), keymetric(19)]);%vscoop (no pt 4 here)
% pairedmat(4,2) = nan;%vsdefect - none such trials 
%pt6 below
pairedmat(4,1) = mean([keymetric(20-skipPT05), keymetric(21-skipPT05)]);%vscoop 
pairedmat(4,2) = mean([keymetric(22-skipPT05), keymetric(23-skipPT05)]);%vsdefect 
%AT adding below 2/11/22 to include data from pt4
if skippt05 == 0
    pairedmat(5,1) = mean([keymetric(15-skipPT05), keymetric(16-skipPT05)]);%vscoop 
    pairedmat(5,2) = keymetric(17-skipPT05);%vsdefect 
elseif skippt05 == 0.5
    pairedmat(5,1) = keymetric(16-skipPT05);%vscoop 
    pairedmat(5,2) = keymetric(17-skipPT05);%vsdefect 
end

input1 = pairedmat(:,1); input2 = pairedmat(:,2);
pairedplotting_helperfx(input1, input2, xlabel1text, xlabel2text, ylabeltext, Ylim)






% below is plot comparing cooperates to defects for whatever key variable
index = find(opponent=="human");
input1 = keymetric(index);
index = find(opponent=="computer");
input2 = keymetric(index);
xlabel1text = 'vsHuman'; 
xlabel2text = 'vsComp';
ylabeltext = ylabeltext; %, 'FontSize', 14);
notpairedplotting_helperfx(input1, input2, xlabel1text, xlabel2text, ylabeltext, Ylim)

%within subject vsperson VS vscomputer
keymetric = keymetric; clear pairedmat
xlabel1text = 'vsHuman'; xlabel2text = 'vsComp';
ylabeltext = ylabeltext; %, 'FontSize', 14);
%pt1 below
pairedmat(1,1) = mean([keymetric(2), keymetric(4), keymetric(5)]);%vsperson 
pairedmat(1,2) = mean([keymetric(1), keymetric(3), keymetric(6)]);%vscomp
%pt2 below
pairedmat(2,1) = mean([keymetric(8), keymetric(10)]);%vsperson
pairedmat(2,2) = mean([keymetric(7), keymetric(9)]);%vscomp
%pt3 below
pairedmat(3,1) = mean([keymetric(12), keymetric(14)]);%vsperson 
pairedmat(3,2) = mean([keymetric(11), keymetric(13)]);%vscomp
%pt5 below
pairedmat(4,1) = keymetric(18-skipPT05);%vsperson (remember, no pt 4 here)
pairedmat(4,2) = keymetric(19-skipPT05);%vscomputer (remember, no pt 4 here)
%pt6 below
pairedmat(5,1) = mean([keymetric(20-skipPT05), keymetric(22-skipPT05)]);%vsperson (remember, no pt 4 here)
pairedmat(5,2) = mean([keymetric(21-skipPT05), keymetric(23-skipPT05)]);%vscomputer (remember, no pt 4 here)
%AT adding below 2/11/22 to include data from pt4
if skippt05 == 0
    pairedmat(6,1) = keymetric(16-skipPT05);%vshuman 
    pairedmat(6,2) = mean([keymetric(15-skipPT05), keymetric(17-skipPT05)]);%vscomputer 
elseif skippt05 == 0.5
    pairedmat(6,1) = keymetric(16-skipPT05);%vshuman 
    pairedmat(6,2) = keymetric(17-skipPT05);%vscomputer 
end



input1 = pairedmat(:,1); input2 = pairedmat(:,2);
pairedplotting_helperfx(input1, input2, xlabel1text, xlabel2text, ylabeltext, Ylim)



2+2

%% more nuanced analyses:
%calculate likelihood of retaliation, session by session

%pseudocode:
%set up some if/then statement to isolate:
%1) subject response when on previous trial, opp defected
%2) subject response when on previous trial, opp cooperated





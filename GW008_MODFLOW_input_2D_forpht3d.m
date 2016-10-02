% MODFLOW_input_print_test_quickstart.m
%
% v3 (10/28/15): notes added to share script

% ############################################################
% Search for: "CUSTOMIZE TO YOUR COMPUTER!!"
% ############################################################


% Prints ASCII files for MODFLOW ba6, dis, and top of lpf input files.
% ===================
% Assumes following:
% ===================
% For ba6 file:
%  - 1D vertical profile domain
%  - Top boundary cond constant head
%  - Bottom boundary cond constant head
%  - Initial head: solve for hyd gradient consistent with Q, K, and h_Bot
%  - constant DZ, DELR, DELC discretization; 
% For dis file:
%  - NPER = 1; ITMUNI = 4 (day); LENUNI = 2 (m)
%  - LAYCBD = 0 for all layers (i.e. no quasi-3D confining bed under layer)
% Hard-coded file names in .nam file:
%  - test.lst     (MODFLOW output file, summary output)
%  - test.pcg (MODFLOW input file for solver)
%  - test.oc  (MODFLOW input file for output specifications)
%  - test.bud     (MODFLOW output file)
%  - testhead.dat (MODFLOW output file)
%  - ibound.dat (MODFLOW output file)
%%%  - testrech.bud (MODFLOW output file)

% (Note: For arrays, must have values for all columns on same line!  This is not req'd for PHT3D arrays....)

clear all, close all; fclose all;

% **************** CUSTOMIZE TO YOUR COMPUTER!! ***************************
% -- MOD_dir: Directory where you keep your MODFLOW test directories
% (you can ignore "isunix" part, that's for my computers)

fl_gcng = 0;  % 1 for gcng

if fl_gcng
    MOD_dir = '/home/gcng/workspace/ModelRuns_scratch/MODFLOW_projects/Minntac/test3';

    % -- MOD_exe: MODFLOW executable 
    % (not actually used, but script will print command for you for your
    % convenience)
    MOD_exe = '/home/gcng/workspace/Models/MODFLOW/Unix/exe/mf2005';

    % *** NEW
    % only use ONE of these lines for slashstr (comment out the other)
    slashstr = '/';  % use this if on Linux
    % slashstr = '\';  % use this if on Windows        
    
else % Patrick
    MOD_dir = 'C:/Hydro_Modeling/MINNTAC_MATLAB_FILES';

    % -- MOD_exe: MODFLOW executable 
    % (not actually used, but script will print command for you for your
    % convenience)
    MOD_exe = 'C:/Hydro_Modeling/mf2005.1_11/bin/mf2005.exe';

    % *** NEW
    % only use ONE of these lines for slashstr (comment out the other)
    slashstr = '/';  % use this if on Linux
    % slashstr = '\';  % use this if on Windows
end


% ************* (end to CUSTOMIZE TO YOUR COMPUTER!!) *********************

% - Domain parameters (for ba6, dis)
nlay = 75;
ncol = 150;   % 

y_scale = 200/nlay; %ratio set by initial harcoded discretization of 200 rows by 400 columns
x_scale = 400/ncol; %ratio set by initial harcoded discretization of 200 rows by 400 columns

domain_len = 220; % meters
domain_bot_elev = -15.45; % meters
domain_top_elev = 0; % top of domain must be at least this elev (include extra space for WT mov't)

% - head boundary conditions
TopHead = [0:(-9.1)/(ncol-1):-9.1];  % head at top boundary (elev nominal at stream bottom)

% - K array (assume isotropic, but can be heterogeneous)
hydcond = ones(nlay,ncol) * 0.56;  % Avg from MW12 S/I/D m/d
hydcond(1:round(81/y_scale),1:round(128/x_scale)) = 0.672; % taken from k_values_Erik_Smith.jpg in google drive
hydcond(1:round(120/y_scale),round(128/x_scale):round(206/x_scale)) = 0.033; % taken from k_values_Erik_Smith.jpg in google drive
hydcond(1:round(81/y_scale),round(207/x_scale):round(296/x_scale)) = 0.672; % taken from k_values_Erik_Smith.jpg in google drive

fl_recharge = 1;  %1: use recharge
hiRate = 0.00084; % m/d determined by Travis' Hydrus model Core 3Dup cummulative bottom flux
loRate = 0.00056;
% -- name of directory with MODFLOW test (create input files in this
% directory)
% *** WARNING!  IT WILL OVERWRITE EXISTING FILES!!! ***
MODtest_dir = 'GW008_MODFLOW_dir';

% -- names of input files (created with this script)
fil_ba6 = 'test_1D.ba6';
fil_dis = 'test_1D.dis';
fil_lpf = 'test_1D.lpf';
fil_nam = 'test.nam';
fil_lmt = 'test.lmt';

% the following are not changed
fil_pcg = 'test.pcg';
fil_oc = 'test.oc';

% the following is created only if recharge stress package is used
if fl_recharge
    fil_rch = 'test.rch';
end    

% - time info
nper = 1;  % number of stress periods
nstp = ones(nper,1); % number of time steps in a stress period
perlen = ones(nper,1)*10; % length of time in a stress period

SsTr_flag = repmat('SS',nper,1);  % 'SS' for steady-state, 'Tr' for transient
if strcmp(SsTr_flag(1,:), 'Tr') % transient
    Ss = 3e-1; % specific storage [1/m] (confined layers, those below uconfined layers)
%     Sy = 0.2;  % specific yield (unconfined layer)
    Sy = Ss * 0.01;  % specific yield (unconfined layer)
    
    nstp(:) = 10;
end
%% =======================================================================

% MODtest_dir0 = fullfile(MOD_dir, MODtest_dir);
MODtest_dir0 = [MOD_dir, slashstr, MODtest_dir];

if ~exist(MODtest_dir0, 'dir')
    mkdir(MOD_dir, MODtest_dir);
end

% LPF parameters
flow_filunit = 34;  % 0: no cell-by-cell flow terms written, 
                    % >0: cell-by-cell flow terms written to this unit number when "SAVE BUDGET", 
                    % <0: cell-by-cell flow for constant-head cells in listing file when "SAVE BUDGET"
                    
% - discretization

DELR = domain_len / ncol;  DELC = DELR;
dz = (domain_top_elev - domain_bot_elev)/nlay;

IBOUND = ones(nlay, ncol);  % 1: variable head 
IBOUND(1:(200/y_scale),1) = -1 %Constant Head at inward boundary
IBOUND((150/y_scale):(200/y_scale),(400/x_scale)) = -1 %Constant Head at GW008
IBOUND((200/y_scale),1:(400/x_scale)) = 0; %No flow base
IBOUND(1,1:(400/x_scale)) = 1
IBOUND(round(199/y_scale),1:round(390/x_scale)) = 0; %no flow base
IBOUND(round(198/y_scale),1:round(380/x_scale)) = 0; %no flow base
IBOUND(round(197/y_scale),1:round(370/x_scale)) = 0; %no flow base
IBOUND(round(196/y_scale),1:round(360/x_scale)) = 0; %no flow base
IBOUND(round(195/y_scale),1:round(350/x_scale)) = 0; %no flow base
IBOUND(round(194/y_scale),1:round(340/x_scale)) = 0; %no flow base
IBOUND(round(193/y_scale),1:round(330/x_scale)) = 0; %no flow base
IBOUND(round(192/y_scale),1:round(320/x_scale)) = 0; %no flow base
IBOUND(round(191/y_scale),1:round(310/x_scale)) = 0; %no flow base
IBOUND(round(190/y_scale),1:round(300/x_scale)) = 0; %no flow base
IBOUND(round(189/y_scale),1:round(290/x_scale)) = 0; %no flow base
IBOUND(round(188/y_scale),1:round(280/x_scale)) = 0; %no flow base
IBOUND(round(187/y_scale),1:round(270/x_scale)) = 0; %no flow base
IBOUND(round(186/y_scale),1:round(260/x_scale)) = 0; %no flow base
IBOUND(round(185/y_scale),1:round(250/x_scale)) = 0; %no flow base
IBOUND(round(184/y_scale),1:round(240/x_scale)) = 0; %no flow base
IBOUND(round(183/y_scale),1:round(230/x_scale)) = 0; %no flow base
IBOUND(round(182/y_scale),1:round(220/x_scale)) = 0; %no flow base
IBOUND(round(181/y_scale),1:round(210/x_scale)) = 0; %no flow base
IBOUND(round(180/y_scale),1:round(200/x_scale)) = 0; %no flow base
IBOUND(round(179/y_scale),1:round(190/x_scale)) = 0; %no flow base
IBOUND(round(178/y_scale),1:round(180/x_scale)) = 0; %no flow base
IBOUND(round(177/y_scale),1:round(170/x_scale)) = 0; %no flow base
IBOUND(round(176/y_scale),1:round(160/x_scale)) = 0; %no flow base
IBOUND(round(175/y_scale),1:round(150/x_scale)) = 0; %no flow base
IBOUND(round(174/y_scale),1:round(140/x_scale)) = 0; %no flow base
IBOUND(round(173/y_scale),1:round(130/x_scale)) = 0; %no flow base
IBOUND(round(172/y_scale),1:round(120/x_scale)) = 0; %no flow base
IBOUND(round(171/y_scale),1:round(110/x_scale)) = 0; %no flow base
IBOUND(round(170/y_scale),1:round(100/x_scale)) = 0; %no flow base
IBOUND(round(169/y_scale),1:round(90/x_scale)) = 0; %no flow base
IBOUND(round(168/y_scale),1:round(80/x_scale)) = 0; %no flow base
IBOUND(round(167/y_scale),1:round(70/x_scale)) = 0; %no flow base
IBOUND(round(166/y_scale),1:round(60/x_scale)) = 0; %no flow base
IBOUND(round(165/y_scale),1:round(50/x_scale)) = 0; %no flow base
IBOUND(round(164/y_scale),1:round(40/x_scale)) = 0; %no flow base
IBOUND(round(163/y_scale),1:round(30/x_scale)) = 0; %no flow base
IBOUND(round(162/y_scale),1:round(20/x_scale)) = 0; %no flow base
IBOUND(round(161/y_scale),1:round(10/x_scale)) = 0; %no flow base



initHead = repmat(TopHead, nlay, 1);

% if fl_recharge == 1 && ~isempty(find(IBOUND(1,:) < 0,1))
%    fprintf('Error! Top layer must be variable head (IBOUND > 0) to receive recharge! Exiting...\n');
%    return
% end


% get layer elevations
botm = [domain_bot_elev: dz: domain_top_elev]';
botm = botm(end-1: -1: 1);
top_sys = domain_top_elev;  % if head > top, this is confined condition


% -- pcg and oc files are not changed with this script
% fil_pcg_0 = fullfile(MODtest_dir0, fil_pcg);
fil_pcg_0 = [MODtest_dir0, slashstr, fil_pcg];
fid = fopen(fil_pcg_0, 'wt');
fprintf(fid, '# Preconditioned conjugate-gradient package\n');
fprintf(fid, '        50        30         1      MXITER, ITER1, NPCOND\n');
fprintf(fid, '  0000.001      .001        1.         2         1         1      1.00\n');
fprintf(fid, ' HCLOSE,      RCLOSE,    RELAX,    NBPOL,     IPRPCG,   MUTPCG    damp\n');
fclose(fid);

% fil_oc_0 = (MODtest_dir0, fil_oc);
fil_oc_0 = [MODtest_dir0, slashstr, fil_oc];
fid = fopen(fil_oc_0, 'wt');
fprintf(fid, 'HEAD PRINT FORMAT 20\n');
fprintf(fid, 'HEAD SAVE UNIT 51\n');
fprintf(fid, 'COMPACT BUDGET AUX\n');
fprintf(fid, 'IBOUND SAVE UNIT 52\n');
for per_i = 1:nper
    for stp_i = 1:nstp(per_i)
        fprintf(fid, 'PERIOD %d STEP %d\n', per_i, stp_i);
        if stp_i == nstp
            fprintf(fid, '   PRINT HEAD\n');
            fprintf(fid, '   SAVE IBOUND\n');
            fprintf(fid, '   PRINT BUDGET\n');
        end
        fprintf(fid, '   SAVE HEAD\n');
        fprintf(fid, '   SAVE BUDGET\n');
    end
end
fclose(fid);

% -- ba6 file
% fil_ba6_0 = fullfile(MODtest_dir0, fil_ba6);
fil_ba6_0 = [MODtest_dir0, slashstr, fil_ba6];
fmt1 = [repmat('%4d', 1, ncol), '\n']; % for IBOUND (for x-section, ow need nrow x ncol for each layer)
fmt2 = [repmat('%7g ', 1, ncol), '\n']; % for initHead

fid = fopen(fil_ba6_0, 'wt');
fprintf(fid, '# basic package file --- x-section %d layers, %d column\n', nlay, ncol);
fprintf(fid, 'XSECTION\n');
fprintf(fid, 'INTERNAL          1 (FREE)  3\n');
fprintf(fid, fmt1, IBOUND');
fprintf(fid, '    999.99  HNOFLO\n');
fprintf(fid, 'INTERNAL          1.0 (FREE)  3\n');
fprintf(fid, fmt2, initHead');
fclose(fid);

% -- dis file
% for now, assume these inputs:
nrow = 1;  % 1 for x-section
itmuni = 4;  % time units, 4: day
lenuni = 2;  % length units, 2: m
LAYCBD = zeros(1,nlay);  % 0: no quasi-3D confining bed under layer
tsmult = 1;  % multiplier for length of successive time steps

% fil_dis_0 = fullfile(MODtest_dir0, fil_dis);
fil_dis_0 = [MODtest_dir0, slashstr, fil_dis];
fmt3 = [repmat(' %d', 1, nlay), '\n']; % for LAYCBD

fid = fopen(fil_dis_0, 'wt');
fprintf(fid, ' %d  %d  %d  %d  %d  %d ', nlay, nrow, ncol, nper, itmuni, lenuni);
fprintf(fid, '  NLAY,NROW,NCOL,NPER,ITMUNI,LENUNI\n');
fprintf(fid, fmt3, LAYCBD);
fprintf(fid, 'CONSTANT %14g   DELR\n', DELR);
fprintf(fid, 'CONSTANT %14g   DELC\n', DELC);
fprintf(fid, 'CONSTANT %14g   TOP of system\n', top_sys);
for ii = 1: nlay
    fprintf(fid, 'CONSTANT %14g   Layer BOTM layer %4d\n', botm(ii), ii);
end
for ii = 1: nper
    fprintf(fid, ' %g %d %g %c%c        PERLEN, NSTP, TSMULT, Ss/Tr (stress period %4d)\n', ...
        perlen(ii), nstp(ii), tsmult, SsTr_flag(ii,:), ii);
end
fclose(fid);


% -- top of lpf file (use Ksat_realization or other to create Ksat inputs)
% assumed input values
hdry = 1e30;  % head assigned to dry cells
nplpf = 0;    % number of LPF parameters (if >0, key words would follow)
laytyp = ones(nlay,1);  % laytyp(1) = 1;  % top is "covertible", rest are "confined"
layave = zeros(nlay,1);  % harmonic mean for interblock transmissivity
chani = ones(nlay,1);   % constant horiz anisotropy mult factor (for each layer)
layvka = zeros(nlay,1);  % flag 0: vka is vert K; >0 is vertK/horK ratio
laywet = zeros(nlay,1);  % flag 0: wetting off


fmt1 = repmat('%2d', 1, nlay);

% fil_lpf_0 = fullfile(MODtest_dir0, fil_lpf);
fil_lpf_0 = [MODtest_dir0, slashstr, fil_lpf];
fid = fopen(fil_lpf_0, 'wt');
fprintf(fid, '# LPF package inputs\n');
fprintf(fid, '%d %g %d    ILPFCB,HDRY,NPLPF\n', flow_filunit, hdry, nplpf);
fprintf(fid, [fmt1, '     LAYTYP\n'], laytyp);
fprintf(fid, [fmt1, '     LAYAVE\n'], layave);
fprintf(fid, [fmt1, '     CHANI \n'], chani);
fprintf(fid, [fmt1, '     LAYVKA\n'], layvka);
fprintf(fid, [fmt1, '     LAYWET\n'], laywet);


% -- Write HKSAT and Ss, Sy (if Tr) in .lpf file
format0 = [repmat(' %4.2f', 1, ncol), '\n'];
format1 = [repmat(' %4.2e', 1, ncol), '\n'];
% loop thru layers (different entry for each layer)
for lay = 1: nlay
    fprintf(fid, 'INTERNAL   1.000E-00 (FREE)    0            HY layer  %d\n', lay);
    fprintf(fid, format0, hydcond(lay,:));

    fprintf(fid, 'INTERNAL   1.000E-00 (FREE)    0            VKA layer  %d\n', lay);
    fprintf(fid, format0, hydcond(lay,:));
    
    if strcmp(SsTr_flag, 'Tr')
        if laytyp(lay) == 0 % confined
            fprintf(fid, 'INTERNAL   1.000E-00 (FREE)    0            Ss layer  %d\n', lay);
            fprintf(fid, format1, ones(ncol,1)*Ss);
        elseif laytyp(lay) == 1 % convertible, i.e. unconfined
            fprintf(fid, 'INTERNAL   1.000E-00 (FREE)    0            Sy layer  %d\n', lay);
            fprintf(fid, format1, ones(ncol,1)*Sy);
        end
    end
end
fprintf(fid, '\n');
fclose(fid);


% -- .rch file
if fl_recharge
    NRCHOP = 3;  % option 3: apply recharge rate to highest active cell in column
    IRCHBD = -1; % <0 to skip writing cell-by-cell flow terms to file
    INRECH = 1; % >0: read in recharge values below (<0: use recharge from previous stress period)
%     fil_rch_0 = fullfile(MODtest_dir0, fil_rch);
    fil_rch_0 = [MODtest_dir0, slashstr, fil_rch];
    fid = fopen(fil_rch_0, 'wt');
    fprintf(fid, '# recharge package \n');  
    fprintf(fid, '    %2d    %2d        NRCHOP,IRCHBD\n', NRCHOP,IRCHBD); % 
    for per_i = 1: nper
        fprintf(fid, '    %2d              INRECH\n', INRECH); % 
        % fprintf(fid, 'CONSTANT %14g   RECH (PERIOD %d) \n', rch_rate, per_i);
        SpatVarRechRate = ones(ncol,1);
        SpatVarRechRate(1:nearest(296/x_scale)) = hiRate;
        SpatVarRechRate(nearest(296/x_scale):nearest(400/x_scale)) = loRate;
        fprintf(fid, 'INTERNAL   1.0 (FREE) 0         recharge rate  \n');
        fprintf(fid, format1, SpatVarRechRate);
    end
    fclose(fid);
end


% -- .nam file with full paths
% fil_nam_0 = fullfile(MODtest_dir0, fil_nam);
fil_nam_0 = [MODtest_dir0, slashstr, fil_nam];
fid = fopen(fil_nam_0, 'wt');
% fprintf(fid, 'LIST          7 %s \n', fullfile(MODtest_dir0, 'test.lst')); % MODFLOW output file
fprintf(fid, 'LIST          7 %s \n', [MODtest_dir0, slashstr, 'test.lst']); % MODFLOW output file
fprintf(fid, 'BAS6          8 %s \n', fil_ba6_0);
fprintf(fid, 'LPF          11 %s \n', fil_lpf_0);
if fl_recharge
%     fprintf(fid, 'RCH          18 %s \n', fullfile(MODtest_dir0, fil_rch));
    fprintf(fid, 'RCH          18 %s \n', [MODtest_dir0, slashstr, fil_rch]);
end
% fprintf(fid, 'PCG          19 %s \n', fullfile(MODtest_dir0, fil_pcg));
fprintf(fid, 'PCG          19 %s \n', [MODtest_dir0, slashstr, fil_pcg]);
% fprintf(fid, 'OC           22 %s \n', fullfile(MODtest_dir0, fil_oc));
fprintf(fid, 'OC           22 %s \n', [MODtest_dir0, slashstr, fil_oc]);
fprintf(fid, 'DIS          10 %s \n', fil_dis_0);
% fprintf(fid, 'LMT6         66 %s \n', fullfile(MODtest_dir0, fil_lmt));
fprintf(fid, 'LMT6         66 %s \n', [MODtest_dir0, slashstr, fil_lmt]);
fprintf(fid, 'DATA(BINARY) 34 %s \n', fullfile(MODtest_dir0, 'test.bud')); % MODFLOW output file
% if fl_recharge
%     fprintf(fid, 'DATA(BINARY) 35 %s \n', fullfile(MODtest_dir0, 'testrech.bud')); % MODFLOW output file
% end
% fprintf(fid, 'DATA(BINARY) 51 %s \n', fullfile(MODtest_dir0, 'testhead.dat')); % MODFLOW output file
fprintf(fid, 'DATA(BINARY) 51 %s \n', [MODtest_dir0, slashstr, 'testhead.dat']); % MODFLOW output file
% fprintf(fid, 'DATA         52 %s \n', fullfile(MODtest_dir0, 'ibound.dat')); % MODFLOW output file
fprintf(fid, 'DATA         52 %s \n', [MODtest_dir0, slashstr, 'ibound.dat']); % MODFLOW output file

fclose(fid);


% -- .lmt file with full paths (this is the binary flow file needed by
% reactive transport model PHT3D)
% fil_lmt_0 = fullfile(MODtest_dir0, fil_lmt);
fil_lmt_0 = [MODtest_dir0, slashstr, fil_lmt];
fid = fopen(fil_lmt_0, 'wt');
fprintf(fid, '# Link-MT3DMS file for MODFLOW-2000, generated by PMWIN \n');
% fprintf(fid, 'OUTPUT_FILE_NAME   %s\n', fullfile(MODtest_dir0, 'test.flo')); % MODFLOW output file
fprintf(fid, 'OUTPUT_FILE_NAME   %s\n', [MODtest_dir0, slashstr, 'test.flo']); % MODFLOW output file
fprintf(fid, 'OUTPUT_FILE_UNIT   333 \n');
fprintf(fid, 'OUTPUT_FILE_HEADER Standard \n');
fprintf(fid, 'OUTPUT_FILE_FORMAT Unformatted \n');
fclose(fid);


% -- Print instructions to run model:
fprintf('MODFLOW input files created!  \n');
fprintf('To run model: \n');
fprintf('(from %s or whichever directory you want output files in)\n', MODtest_dir0);
fprintf('\nprompt> %s\n\n', [MOD_exe, ' ', fil_nam_0]);


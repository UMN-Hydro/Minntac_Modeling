function [mob_eq_comp, mob_eq_ic_z, mob_eq_extra_z, min_eq_comp, min_eq_ic_z, ...
    catex_comp, catex_ic_z, surf_comp, surf_ic_z, surf_par, surf_cpl, surf_calc_type] = ...
    MW12_Minntac_InitCond_chem_red_wFeS(phrq_sim_dir, phrq_exe, use_file_databas, por, tempC) 

% Uses the following functions:
%   - Alk2DIC() (only needed if data is for alkalinity)
%       Converts alkalinity to total C(4) 
%   - generate_ic_PHREEQC_f_051915a()
%       Generates equilibrated ic using PHREEQC
%   - go_PHREEQC2_3()
%       Called by generate_ic_PHREEQC_f_051915a(), runs PHREEQC

% % - ONLY uncomment below if running as script
 fl_gcng = 0; % 1 for gcng, 0 for Patrick
 if fl_gcng
     sim_dir = 'C:\Hydro_Modeling\pht3d_dir\';
     phrq_exe = 'C:\Hydro_Modeling\phreeqc.exe';
     use_file_databas = 'C:\Hydro_Modeling\pht3d_dir\pht3d_datab2.dat';    
 else
     sim_dir = 'C:\Hydro_Modeling\pht3d_dir\';
     phrq_exe = 'C:\Hydro_Modeling\phreeqc.exe';
     use_file_databas = 'C:\Hydro_Modeling\pht3d_dir\pht3d_datab.dat_160331';
 end
 tempC = 10.0; 
 por = 0.39;
 phrq_sim_dir = sim_dir;


% % **** To run as script: comment out 1st line for function, uncomment below
% % block
% SC_dir = '/home/gcng/workspace/matlab_files/SecondCreek';
% addpath(SC_dir);
% sim_dir = '/home/gcng/workspace/ModelRuns_scratch/PHT3D_projects/SecondCreek/test1_1D/';
% phrq_sim_dir = sim_dir;
% 
% % you may not need this
% matlab_dir = '/home/gcng/workspace/matlab_files/my_toolbox/PHT3D_functions';  
% addpath(matlab_dir);
% % ----------------

% *** CUSTOMIZE: **********************************************************
% This script uses PHREEQC to prepare (equilibrate and charge-balance)
% different aqueous solutions for use in initial and boundary conditions.  
% Different aq solutions could correspond to different parts of the domain,
% e.g. Cell 2 water, bulk domain groundwater, and recharge water.  These
% different aq solutions will be called 'zones'.
%  
% Use this section to specify the number of zones and to comment on what
% each zone represents.
nzones = 4;
% ** comment here to describe 'zones'
% zone 1: cell 2 water
% zone 2: domain ic water
% zone 3: recharge water down-gradient of dike
% zone 4: recharge water thru dike
zone_names = {'Cell', 'domain ic', 'recharge down-grad', 'recharge dike'};
% *************************************************************************

% -- SET MODEL CHEMISTRY: initial trial values
% (Units -- aq: mol/L_w, user-defined immob (e.g. bacteria, napl): mol/L_w, 
% minerals (and gases): mol/L_v, exchangers and surfaces: mol/L_v)

% - mobile equil components
% (ic from PHREEQC batch run)    
n_mob_eq_max = 50;
mob_eq_comp = cell(n_mob_eq_max,1);

mob_eq_obs_z = zeros(nzones, n_mob_eq_max);
convert_obs2model = zeros(1, n_mob_eq_max); % convert obs units to mol/L
mob_eq_extra_z = cell(nzones, n_mob_eq_max);
ii = 0; 
ii=ii+1; mob_eq_comp{ii} = 'C(4)'; % ***********************************
% [alkalinity as g/L CaCO3] -> [alkalinity as mol/L equivalence]: divide by
% 50g/eq (because CaCO3 has 100 g/mol molar mass, 1 mol CO3^{2-} / 2 eq)
Ka1 = 10^-6.3; Ka2 = 10^-10.3; KH = 10^-14;
pH = 6.92;
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = Alk2DIC(0.25/50, pH, Ka1, Ka2, KH); % meq/L (alk = 207.5 mg/L as CaCO3 at Cell 1, 307.5 mg/L as CaCO3 at Cell 2) 
zz = zz + 1; mob_eq_obs_z(zz,ii) = Alk2DIC(0.15/50, pH, Ka1, Ka2, KH); % eq/L (alk = 150 mg/L as CaCO3 at MW12)
zz = zz + 1; mob_eq_obs_z(zz,ii) = Alk2DIC(0.15/50, pH, Ka1, Ka2, KH); % eq/L (alk = 150 mg/L as CaCO3 at MW12)
zz = zz + 1; mob_eq_obs_z(zz,ii) = Alk2DIC(0.15/50, pH, Ka1, Ka2, KH); % eq/L (alk = 150 mg/L as CaCO3 at MW12)
convert_obs2model(ii) = 1;
ii=ii+1; mob_eq_comp{ii} = 'C(-4)'; % ***********************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L 
convert_obs2model(ii) = 1/1e3/16.0;
% ii=ii+1; mob_eq_comp{ii} = 'Ca'; % ***********************************
% zz = 0;
% zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L
% zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L
% zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L
% zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L
% convert_obs2model(ii) = 1/1000/40.08;
ii=ii+1; mob_eq_comp{ii} = 'Cl'; % ***********************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 141.557; % mg/L
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L
convert_obs2model(ii) = 1/1000/35.453;
% ii=ii+1; mob_eq_comp{ii} = 'F'; % ***********************************
% zz = 0;
% zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L
% % mob_eq_extra_z{zz,ii} = 'charge';
% zz = zz + 1; mob_eq_obs_z(zz,ii) = 1; % mg/L
% mob_eq_extra_z{zz,ii} = 'charge';
% zz = zz + 1; mob_eq_obs_z(zz,ii) = 1; % mg/L
% mob_eq_extra_z{zz,ii} = 'charge';
% zz = zz + 1; mob_eq_obs_z(zz,ii) = 1; % mg/L
% mob_eq_extra_z{zz,ii} = 'charge';
% convert_obs2model(ii) = 1/1000/35.453;
ii=ii+1; mob_eq_comp{ii} = 'Fe(2)'; % ***********************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % ug Fe/L
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % ug Fe/L
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % ug Fe/L
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % ug Fe/L
convert_obs2model(ii) = 1/1e6/55.847;
% ii=ii+1; mob_eq_comp{ii} = 'K'; % ***********************************
% zz = 0;
% zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L
% zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L
% zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L
% zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L
% convert_obs2model(ii) = 1/1e3/39.102;
% ii=ii+1; mob_eq_comp{ii} = 'Mg'; % ***********************************
% zz = 0;
% zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L
% zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L
% zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L
% zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L
% convert_obs2model(ii) = 1/1e3/24.312;
% ii=ii+1; mob_eq_comp{ii} = 'Mn(2)'; % ***********************************
% zz = 0;
% zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % ug/L
% zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % ug/L
% zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % ug/L
% zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % ug/L
% convert_obs2model(ii) = 1/1e6/54.938;
ii=ii+1; mob_eq_comp{ii} = 'Na'; % ***********************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 1; % mg/L 'charge' denotes will balance with anions
mob_eq_extra_z{zz,ii} = 'charge';
zz = zz + 1; mob_eq_obs_z(zz,ii) = 1; % mg/L
mob_eq_extra_z{zz,ii} = 'charge';
zz = zz + 1; mob_eq_obs_z(zz,ii) = 1; % mg/L
mob_eq_extra_z{zz,ii} = 'charge';
zz = zz + 1; mob_eq_obs_z(zz,ii) = 1; % mg/L
mob_eq_extra_z{zz,ii} = 'charge';
convert_obs2model(ii) = 1/1e3/22.9898;
ii=ii+1; mob_eq_comp{ii} = 'O(0)'; % ***********************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 9; % mg/L DO (saturation is 9mg/L at 20C) 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L DO  
zz = zz + 1; mob_eq_obs_z(zz,ii) = 9; % mg/L DO  
zz = zz + 1; mob_eq_obs_z(zz,ii) = 9; % mg/L DO  
convert_obs2model(ii) = 1/1e3/16.0;
ii=ii+1; mob_eq_comp{ii} = 'S(-2)'; % ***********************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg S/L
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg S/L
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg S/L
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg S/L
convert_obs2model(ii) = 1/1e3/32.064;
ii=ii+1; mob_eq_comp{ii} = 'S(6)'; % ***********************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 935.67; % mg SO4/L
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg SO4/L
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg SO4/L
% thru tailings: 15 metric tons / mi2 / wk (= 0.0726 mg/L)
Rech = 0.0114; % recharge [m/d]
L = 15/7/1609.34^2*1e9; % loading [mg/m2/d]
conc = L/Rech/1e3; % mg S/L
zz = zz + 1; mob_eq_obs_z(zz,ii) = conc; % mg S/L
convert_obs2model(ii) = 1/1e3/(32.064+4*16.0);
% ii=ii+1; mob_eq_comp{ii} = 'Si'; % ***********************************
% zz = 0;
% zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L
% zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L
% zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L
% zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L
% convert_obs2model(ii) = 1/1e3/28.0843;
ii=ii+1; mob_eq_comp{ii} = 'pH'; % ***********************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 6.92; % 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 6.92; % 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 6.92; % 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 6.92; % 
convert_obs2model(ii) = 1;
ii=ii+1; mob_eq_comp{ii} = 'pe'; % ***********************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 14; % 
zz = zz + 1; mob_eq_obs_z(zz,ii) = -7; % 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 14; % 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 14; % 
convert_obs2model(ii) = 1;
n_mob_eq = ii;
mob_eq_comp = mob_eq_comp(1:n_mob_eq);
mob_eq_obs_z = mob_eq_obs_z(:,1:n_mob_eq);
mob_eq_extra_z = mob_eq_extra_z(:,1:n_mob_eq);
convert_obs2model = convert_obs2model(1:n_mob_eq);
mob_eq_obs_z = mob_eq_obs_z .* repmat(convert_obs2model, nzones, 1);


% - mineral eq components
% (only allowed for porewater solution)
n_min_eq_max = 10;
min_eq_comp = cell(n_min_eq_max,1);
min_eq_ic_z = zeros(nzones,n_min_eq_max); 
ii = 0;
% ii=ii+1; min_eq_comp{ii} = 'Fe(OH)3(a)'; % ***********************************
% min_eq_ic_z(zz,ii) = 10;  % arbitrary conc
% ii=ii+1; min_eq_comp{ii} = 'Gibbsite'; % ***********************************
% min_eq_ic_z(zz,ii) = 30e-6*rho_b; % (mol/Lv)
% ii=ii+1; min_eq_comp{ii} = 'Goethite'; % ***********************************
% min_eq_ic_z(zz,ii) = 50e-6*rho_b; % (mol/Lv)
% ii=ii+1; min_eq_comp{ii} = 'Pyrolusite'; % ***********************************
% min_eq_ic_z(zz,ii) = .25e-6*rho_b;  % 
ii=ii+1; min_eq_comp{ii} = 'Siderite'; % ***********************************
zz = 0;
zz = zz + 1; min_eq_ic_z(zz,ii) = 0;  % 0 mineral in Cell 2 water
zz = zz + 1; min_eq_ic_z(zz,ii) = 10;  % arbitrary conc
zz = zz + 1; min_eq_ic_z(zz,ii) = 0;  % arbitrary conc
ii=ii+1; min_eq_comp{ii} = 'FeS(ppt)'; % ***********************************
zz = 0;
zz = zz + 1; min_eq_ic_z(zz,ii) = 0;  % 0 mineral in Cell 2 water
zz = zz + 1; min_eq_ic_z(zz,ii) = 10;  % arbitrary conc
zz = zz + 1; min_eq_ic_z(zz,ii) = 0;  % arbitrary conc
% ii=ii+1; min_eq_comp{ii} = 'Rhodochrosite'; % ***********************************
% ii=ii+1; min_eq_comp{ii} = 'Vivianite'; % ***********************************
n_min_eq = ii;
min_eq_comp = min_eq_comp(1:n_min_eq);
min_eq_ic_z = min_eq_ic_z(:,1:n_min_eq); 


% - catex components (catex_exch is exchanger site)
n_catex_exch_max = 10;
catex_exch_comp = cell(n_catex_exch_max,1);
catex_exch_ic_z = zeros(nzones, n_catex_exch_max); 
ii = 0;
n_catex_exch = ii;
catex_exch_comp = catex_exch_comp(1:n_catex_exch);
catex_exch_ic_z = catex_exch_ic_z(:,1:n_catex_exch); 

% - catex components (catex is exchanger sorbed compound)
% (ic from PHREEQC batch run)
n_catex_max = 10;
catex_comp = cell(n_catex_max,1);
ii = 0;
n_catex = ii;
catex_comp = catex_comp(1:n_catex);


% - surface complexation components
% **** WARNING: not set up for surfaces w/
%               multiple sites (would need to enter only 1 surf_calc_type
%               per surface in _ph file)
n_surf_max = 10;
surf_comp = cell(n_surf_max,1);
% surf_ic: number of sites [mol/Lv] if not coupled to phase or reactant,
%          (number of sites per mol x porosity) [mol/mol] if coupled
surf_ic_z = zeros(nzones,n_surf_max); 
% surf_par(1,ii): SurfArea ([m2/g] if not coupled to phase or reactant,
%              [mol/m2] if coupled to phase or reactant)
% surf_par(2,ii): mass ([g/L] if not coupled to phase or reactant,
%              include but ignored if coupled to phase or reactant)
surf_par = zeros(2,n_surf_max); 
% surf_cpl{1,ii}: Empty if not coupled to phase or reactant
% surf_cpl{1,ii}: Name of pure phase or kin reactant coupled to
% surf_cpl{2,ii}: 'equilibrium_phase' to couple to pure phase, or
%                    'kinetic_reactant' to couple to kinetic reactant
surf_cpl = cell(2,n_surf_max); 
surf_calc_type = ''; 
ii = 0;
n_surf = ii;
surf_comp = surf_comp(1:n_surf);
surf_ic_z = surf_ic_z(:,1:n_surf); 
surf_par = surf_par(:,1:n_surf); 

% (done setting model chemical components)

% - select output list
phrq_sel_outlist.total = [mob_eq_comp(~strncmp('p',mob_eq_comp,1))];
phrq_sel_outlist.mol = catex_comp;
phrq_sel_outlist.equilphase = min_eq_comp;
phrq_sel_outlist.si = []; phrq_sel_outlist.gas = [];    


% -- Get initial conditions w/ PHREEQC
%     % (re-)run PHREEQC to equilibrate ic conc for each zone
%     fprintf('Have not finished code of fl_init_phrq=1!  Need to do for each zone!  Exiting...\n');
%     return
mob_eq.name = mob_eq_comp; mob_eq.extra_z = mob_eq_extra_z;


if n_catex == 0 % not cation exchange components
    catex = []; 
else
    % assume exchanger is same for all zones
    catex.name = catex_exch_comp; 
    catex.cec = catex_exch_ic / por;
end
fl_force_redox = 1;  % force redox equilibration, not just charge balance


% -- Set file names
phrq_infil = fullfile(phrq_sim_dir, 'batch_phrq.dat'); 
phrq_outfil = fullfile(phrq_sim_dir, 'batch_phrq.out'); 
phrq_tblfil = fullfile(phrq_sim_dir, 'batch_phrq.tbl'); 
phrq_sel_outfil = fullfile(phrq_sim_dir, 'batch_phrq.sel'); 

% set database file (must be in sim_dir with file name 'pht3d_datab.dat')
file_databas = fullfile(phrq_sim_dir, 'pht3d_datab.dat');
copyfile(use_file_databas, file_databas);


% -- Check to make sure you don't over-write any output files
% - go to simulation dir (check to make sure not over-writing, no current programs running)
if ~exist(phrq_sim_dir, 'dir')
    mkdir(phrq_sim_dir);
end
%if exist([phrq_sel_outfil, '1'], 'file')    
%    fprintf('phrq_sim_dir has .sel file(s)!  Could be job running or unsaved job there.  Exiting...\n');
%    fprintf('(phrq_sim_dir %s) \n', phrq_sim_dir);
%    mob_eq_comp = []; mob_eq_ic_z = []; mob_eq_extra_z = []; min_eq_comp = []; min_eq_ic_z= [];
%    catex_comp = [];catex_ic_z = [];surf_comp = [];surf_ic_z = [];surf_par = [];surf_cpl = [];surf_calc_type = [];    
%    return
%end

% -- Run PHREEQC to equilibrate and charge-balance aq solutions for each zone
units = 'mol/kgw';
mob_eq_ic_z = zeros(nzones, n_mob_eq);
catex_ic_z = nan(nzones, n_catex);
for zz = 1:nzones  % loop through different zones
    fprintf('\n--- Zone %d (%s) ---', zz, zone_names{zz});
    mob_eq.ic = mob_eq_obs_z(zz,:);
    mob_eq.extra = mob_eq_extra_z(zz,:);
    phase.name = min_eq_comp; phase.ic = min_eq_ic_z(zz,:);
    phase.si = zeros(n_min_eq,1); 
    [mob_eq_ic_z(zz,:), b] = generate_ic_PHREEQC_f_062215a_2(...
        phrq_exe, phrq_infil, phrq_outfil, file_databas, phrq_sel_outfil, phrq_sel_outlist, ...
        por, tempC, units, ...
        mob_eq, phase, catex, ...
        mob_eq_comp, catex_comp, fl_force_redox, ...
        phrq_tblfil);
    if ~isempty(b)
        catex_ic_z(zz,:) = b;
    end
    
    movefile(phrq_infil, [phrq_infil, num2str(zz)]);
    movefile(phrq_outfil, [phrq_outfil, num2str(zz)]);
    movefile(phrq_sel_outfil, [phrq_sel_outfil, num2str(zz)]);
    movefile(phrq_tblfil, [phrq_tblfil, num2str(zz)]);
    
end



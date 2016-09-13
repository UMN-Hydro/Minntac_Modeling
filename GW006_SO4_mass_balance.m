% SO4_mass_balance.m
% 9/12/16
%
% Developed for Minntac project.  Summarizes mass balance results for:
%   - SO4 in through upgradient constant head, constant conc BC
%   - SO4 in through dike recharge
%   - SO4 out through downgradient constant head BC
%   - SO4 loss to Sulfide
%   - Sulfide loss to FeS(ppt)
%
% Uses function read_pht3dout_massbal3 (originally created for Bemidji).
% See read_pht3dout_massbal3.m script for info on output.

clear all, close all

fl_gcng = 0;

if fl_gcng
    % PHT3D_oudir: directory with pht3d.out, pht3d.m, and .MAS files
    % (160404a: MW12 simulation for SME2016)
    PHT3D_outdir = '/home/gcng/workspace/ModelRuns_scratch/PHT3D_projects/Minntac/test1/160404a/';
    slashstr = '/';
else
    % PHT3D_oudir: directory with pht3d.out, pht3d.m, and .MAS files
    PHT3D_outdir = 'C:\Hydro_Modeling\GW006_pht3d_dir'; % enter here
    slashstr = '\';    
end

%% ------------------------------------------------------------------------
if ~strcmp(PHT3D_outdir(end), slashstr)
    PHT3D_outdir = [PHT3D_outdir, slashstr];
end
pht3dout_fil = [PHT3D_outdir, 'pht3d.out'];
pht3dm_fil = [PHT3D_outdir, 'pht3d.m'];
MAS_dir = PHT3D_outdir;

MBpht3d = read_pht3dout_massbal3(pht3dout_fil, pht3dm_fil, MAS_dir);

SO4_ind = find(strcmp(MBpht3d.comp_name, 'S_6'));
if isempty(SO4_ind)
    fprintf('Error! No SO4 output file in directory! Exiting... \n');
    return
end
Sulf_ind = find(strcmp(MBpht3d.comp_name, 'S__2'));
if isempty(Sulf_ind)
    fprintf('Error! No Sulfide output file in directory! Exiting... \n');
    return
end
FeS_ind = find(strcmp(MBpht3d.comp_name, 'FeS_ppt'));
if isempty(FeS_ind)
    fprintf('Error! No FeS(ppt) output file in directory! Exiting... \n');
    return
end

times = MBpht3d.time; % [ntimex1], days

% mass balance: >0 for inputs to soln, <0 for outputs from soln [kmol]
% ([ntime x 1])
SO4_Upgrad = MBpht3d.MassBalIn(:, 1, SO4_ind) + MBpht3d.MassBalOut(:, 1, SO4_ind); % net into soln, should be >0
SO4_Rech = MBpht3d.MassBalIn(:, 3, SO4_ind) + MBpht3d.MassBalOut(:, 3, SO4_ind); % net into soln, should be >0 
SO4_Downgrad = MBpht3d.MassBalIn(:, 2, SO4_ind) + MBpht3d.MassBalOut(:, 2, SO4_ind); % net into soln, should be <0
SO4_Red = MBpht3d.MassBalIn(:, 4, SO4_ind) + MBpht3d.MassBalOut(:, 4, SO4_ind); % net into soln, should be <0

FeSppt = MBpht3d.MassBalIn(:, 4, FeS_ind) + MBpht3d.MassBalOut(:, 4, FeS_ind); % net into soln, should be >0

timesYr = times/365;

h1 = nan(5,1); h2 = nan(5,1); leg = cell(5,1);
SRates_all = nan(length(timesYr)-1, 5);
for ii = 1:5
    switch ii
        case 1
            Y = SO4_Upgrad;
            ti = 'Upgrad In';
            LS = '-'; cv = [0 0 0];
        case 2
            Y = SO4_Rech;
            ti = 'Recharge In';
            LS = '--'; cv = [0 0 0];
        case 3
            Y = -SO4_Downgrad;
            ti = 'Downgrad Out';
            LS = '-.'; cv = [0 0 0];
        case 4
            Y = -SO4_Red;
            ti = 'Reduction Out';
            LS = ':'; cv = [0 0 0];
        case 5
            Y = FeSppt;
            ti = 'Ppt Out';
            LS = ':'; cv = [1 0 0];
    end                   

    subplot(2,1,1)
    h1(ii) = line(timesYr, Y, 'LineStyle', LS, 'Color', cv); hold on,
    leg{ii} = ti;
    xlabel('time (yrs)');
    ylabel('kmol S');
    title('Cumulative sulfate mass balance');

    subplot(2,1,2)
    Y2 = (Y(2:end)-Y(1:end-1))./(timesYr(2:end)-timesYr(1:end-1));
    h2(ii) = line(timesYr(2:end), Y2, 'LineStyle', LS, 'Color', cv); hold on,
    leg{ii} = ti;
    
    SRates_all(:,ii) = Y2;
end

subplot(2,1,1)
legend(h1, leg, 'Location', 'NorthWest');
xlabel('time (yrs)');
ylabel('kmol S');
title('Cumulative sulfate mass');
box on

subplot(2,1,2)
xlabel('time (yrs)');
ylabel('kmol S / yr');
title('Sulfate mass rate');
box on


% -- Summary mass balance over entire simulation period)
fprintf('End time (%5.2f yr) sulfate rates [kmol S/yr]:\n', timesYr(end));
for ii = 1:5
    fprintf(' %18s: %g \n', leg{ii}, SRates_all(end,ii));
end




function MBpht3d = read_pht3dout_massbal3(pht3dout_fil, pht3dm_fil, MAS_dir)
%
% 5/31/14
% Reads in detailed cumulative mass balance results in pht3d.out file; also
% uses pht3d.m to get component names.
%
% 11/28/14: Strangely, .MAS sometimes has many more time printouts.  This
% scripts reads all in, but only extracts thoses corresponding to
% PHT3D_OUTPUT_TIMES.ACN  
%
% 9/12/16: Edits to comments
%
% (note: this mass balance in more detailed than the one in .MAS files)
%
% Output: MBpht3d.
%   BalTitle = {'CONSTANT CONCENTRATION'; 'CONSTANT HEAD'; 'RECHARGE';
%                 'REACTION (PHREEQC)'; 'MASS STORAGE (SOLUTE)'}
%   comp_name = {ncomp,1}
%   time = [ntime x 1]
%   MassBalIn = [ntime,5,ncomp]   (>0, read from pht3d.out)
%   MassBalOut = [ntime,5,ncomp]  (<0, read from pht3d.out)
%   MassAquif = [ntime,ncomp]     (>0, read from *.MAS files) 
%
% Notes:
%   - All results in kmol, but solids are further (incorrectly) multiplied
%     by por_n. kmol units was verified for Bemidji and Minntac runs
%   - Note that files are INCORRECTLY labelled "kg"
%   - All "MassBalOut" are <0
%   - BalTitle: 
%       CONSTANT CONCENTRATION: mass change at constant conc boundary
%           (left-hand boundary)
%       CONSTANT HEAD: mass change at constant head boundary that is not a
%          constant conc boundary (right-hand boundary)
%       RECHARGE: mass change with recharge
%       REACTION (PHREEQC): mass change due to reactions, including both
%          equilibrium (e.g. mineral phases) and kinetic
%       MASS STORAGE (SOLUTE): IN is removal from aquifer storage into the 
%          reaction solution (e.g, for reactions and transport out); OUT is
%          solutes from the reaction solution going into aquifer storage. 
%          Change in MassAquif is: -(OUT MASS STORAGE + IN MASS STORAGE).  
%          (Added because OUT is always <0)
%
% Comments to clarify mass balance output:
%   -- MassAquif ("TOTAL MASS IN AQ") is total mass (kmol) in domain
%   -- MassAquif = -(OUT MASS STORAGE + IN MASS STORAGE) sometimes, but not
%         always.  Thus, use "TOTAL MASS IN AQ" in .MAS file!!
%   -- Also, SHOULD have: 
%           MassAquif(t) = MassAquif(t-1) + ConstConc_NetIn + 
%                           ConstHead_NetIn + Rech_NetIn + Rxn_NetIn
%   -- "DISCREPANCY(%)" in .MAS: Corresponds to pht3d.out "NET" and
%         "DISCREPANCY" results, difference btwn calculated mass change in
%         aquifer and net sinks and sources (sinks & sources refer to
%         transport and rxn)
%   -- .MAS file:
%       -- Mass balance numbers are for a rxn solution, NOT the actual
%          aquifer.  Mass Storage in aquifer is a different entity than the
%          reaction solution.  
%       -- "TOTAL IN" = sum all IN's in pht3d.out (include MASS STORAGE)
%       -- "TOTAL OUT" = sum all OUT's in pht3d.out (include MASS STORAGE)
%       -- "SOURCES" = ConstConc_In + ConstHead_In + Rech_In + Rxn_In
%               (i.e., does not include inputs from Storage)
%       -- "SINKS" = ConstConc_Out + ConstHead_Out + Rech_Out + Rxn_Out 
%               (i.e., does not include output to Storage)
%       -- "TOTAL MASS IN AQ" approx MassAquif(t-1) + ConstConc_NetIn + 
%                           ConstHead_NetIn + Rech_NetIn + Rxn_NetIn
%           (still confuses me why not exactly)
%
%   -- pht3d.out file
%       -- Mass balance numbers are for a rxn solution, NOT the actual
%          aquifer.  Mass Storage in aquifer is a different entity.  
%       -- OUT MASS STORAGE is mass going into storage from rxn solution
%       -- IN MASS STORAGE is mass coming from storage into rxn solution.
%          Note that this can be greater than MassAquif(t-1) because it
%          includes solute sources during that same time period.  Can also
%          be less than MassAquif(t-1), because apparently there can be
%          some storage solute from previous time that is not moved to rxn
%          solution.  CONFUSING!!
%       

% % test as script:
% clear all, close all, fclose all;
% pht3dout_fil = 'C:\cygwin\home\gng\Models\PHT3D\gng_tests_Bemidji\Nvoc_inhib_nalk\053114a\pht3d.out';
% pht3dm_fil = 'C:\cygwin\home\gng\Models\PHT3D\gng_tests_Bemidji\Nvoc_inhib_nalk\053114a\pht3d.m';
% MAS_dir = 'C:\cygwin\home\gng\Models\PHT3D\gng_tests_Bemidji\Nvoc_inhib_nalk\053114a\';
% clear all, close all, fclose all;
% pht3dout_fil = 'C:\cygwin\home\gng\Models\PHT3D\gng_tests_TCE\Soynapl_homog_c3\042214h_p\pht3d.out';
% pht3dm_fil = 'C:\cygwin\home\gng\Models\PHT3D\gng_tests_TCE\Soynapl_homog_c3\042214h_p\pht3d.m';
% MAS_dir = 'C:\cygwin\home\gng\Models\PHT3D\gng_tests_TCE\Soynapl_homog_c3\042214h_p\';
%% ------------------------------------------------------------------------

BalTitle = {'CONSTANT CONCENTRATION'; 'CONSTANT HEAD'; 'RECHARGE'; ...
            'REACTION (PHREEQC)'; 'MASS STORAGE (SOLUTE)'};

% - First get component names
ncomp_max = 50;
comp_name = cell(ncomp_max,1);
fid = fopen(pht3dm_fil, 'r');
ii = 0;
while(1)
    line = fgets(fid);
    if line == -1, break, end
    ii = ii+1;
    ind = find(line == '=');
    comp_name{ii} = line(1:ind-2);
end
fclose(fid);
ncomp = ii;
comp_name = comp_name(1:ncomp);

% - Get cumulative mass balance results
fid = fopen(pht3dout_fil, 'r');
while(1)
    line = fgets(fid);
    ind = regexp(line, 'NUMBER OF TIMES AT WHICH SIMULATION RESULTS ARE SAVED');
    if ~isempty(ind), break, end
end
ind1 = find(line == '='); 
ntime = sscanf(line(ind1+1:end), '%f');
fgets(fid);
time = fscanf(fid, '%f', ntime);

MassBalIn = nan(ntime, length(BalTitle), ncomp);
MassBalOut = nan(ntime, length(BalTitle), ncomp);
for tt = 1: ntime
    for ii = 1: ncomp
        iistr = num2str(ii);
        if ii < 10, iistr = ['0', num2str(ii)]; end
        while(1)
            line = fgets(fid);
            ind = regexp(line, ['FOR COMPONENT NO. ', iistr]);
            if ~isempty(ind), break, end
        end

        % get mass balance
        for jj = 1:18, fgets(fid);  end
        for jj = 1: length(BalTitle)
            line = fgets(fid);
            if isempty(regexp(line, BalTitle{jj}(1:end-1), 'once'))
                fprintf('Error!  Did not find %s in string! Exiting...\n', BalTitle{jj});
                return
            end
            ind1 = find(line == ':'); 
            data = sscanf(line(ind1+1:end), '%f');
            MassBalIn(tt,jj,ii) = data(1);
            MassBalOut(tt,jj,ii) = data(2);
        end
    end
end
fclose(fid);

% - Get total aquifer mass results
MassAquif = nan(ntime, ncomp);
for ii = 1: ncomp
    iistr = num2str(ii);
    if ii < 10, iistr = ['0', num2str(ii)]; end
    
    fid = fopen([MAS_dir, 'PHT3D0', iistr, '.MAS'], 'r');
    for jj =1:2, fgets(fid); end
    data = fscanf(fid, [repmat('%g',1,9), '\n']);
    data = reshape(data,9,length(data)/9);
    fclose(fid);
    MassAquif0 = data(7,:)';
    % make sure only take dates corresp to ACN file
    if ii == 1
        tt2 = zeros(ntime,1);
        for tt = 1: ntime
            b = find(abs(data(1,:)'-time(tt))<1e-3);
            tt2(tt) = b;
        end
        time2 = data(1,tt2)';
    end
    MassAquif(:,ii) = data(7,tt2)';
end


MBpht3d.BalTitle = BalTitle;
MBpht3d.comp_name = comp_name;
MBpht3d.time = time;
MBpht3d.MassBalIn = MassBalIn;
MBpht3d.MassBalOut = MassBalOut;
MBpht3d.MassAquif = MassAquif;


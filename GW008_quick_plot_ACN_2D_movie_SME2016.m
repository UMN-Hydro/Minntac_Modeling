% quickplot
clear all, close all, fclose all;

% **************** CUSTOMIZE TO YOUR COMPUTER!! ***************************
% (you can ignore "isunix" part, that's for my computers)

% - sim_dir: Directory with input files and simulation results
if isunix  
    slashstr = '/';
else
    slashstr = '\';
end
% sim_dir = '/home/gcng/workspace/ModelRuns_scratch/PHT3D_projects/SecondCreek/test2/112915g'; % no slash at end
% sim_dir = '/home/gcng/workspace/ModelRuns_scratch/PHT3D_projects/Minntac/test1'; % no slash at end
% sim_dir = '/home/gcng/workspace/ModelRuns_scratch/PHT3D_projects/Minntac/test1/160310a_DO'; % no slash at end
% sim_dir = '/home/gcng/workspace/ModelRuns_scratch/PHT3D_projects/Minntac/test1/160310b_intenseSO4red'; % no slash at end
% sim_dir = '/home/gcng/workspace/ModelRuns_scratch/PHT3D_projects/Minntac/test1/160331a'; % no slash at end
% sim_dir = '/home/gcng/workspace/ModelRuns_scratch/PHT3D_projects/Minntac/test1/160331f'; % no slash at end
% sim_dir = '/home/gcng/workspace/ModelRuns_scratch/PHT3D_projects/Minntac/test1/160331i'; % no slash at end

% sim_dir = '/home/gcng/workspace/ModelRuns_scratch/PHT3D_projects/Minntac/test1/160401a'; % no slash at end
sim_dir = 'C:\Hydro_Modeling\GW008_pht3d_dir'; % no slash at end
ibound_fil = 'C:\Hydro_Modeling\MINNTAC_MATLAB_FILES\GW008_MODFLOW_dir\ibound.dat';

mytoolboxdir = '/home/gcng/workspace/matlab_files/my_toolbox/';
% ************* (end to CUSTOMIZE TO YOUR COMPUTER!!) *********************

addpath(genpath(mytoolboxdir));

fl_print = 1; % 1: print final figure to .tiff
fl_mov = 1; 
fl_recmov = 1;
fps = 5; % frames per sec (typically 5 is good)
fl_cbar = 1; % include colorbar

fs = 12; % fontsize

domain_len = 224; % meters
domain_bot_elev = -15.45; % meters
domain_top_elev = 0; % top of domain must be at least this elev (include extra space for WT mov't)

% -- what to plot (check pht3d.m in simdir for your options)
n_pcomp_max = 20;
plot_comp = cell(n_pcomp_max,1); 
plot_ti = cell(n_pcomp_max,1); 
xlimv = nan(n_pcomp_max, 2);
sw_clrmp_all = zeros(n_pcomp_max, 1); % colormap code
fl_diff_t1 = zeros(n_pcomp_max,1); 
convert_unit = ones(n_pcomp_max,1); 
unit_ti = repmat({'M'},n_pcomp_max,1);
ii = 0;
ii = ii + 1; plot_comp{ii} = 'Cl'; %xlimv(ii,:) = [0 5];
plot_ti{ii} = 'Chloride'; sw_clrmp_all(ii) = 1;
convert_unit(ii) = 1000*35.453; unit_ti{ii} = 'mg/L'; 
ii = ii + 1; plot_comp{ii} = 'S(6)'; %xlimv(ii,:) = [0 5];
plot_ti{ii} = 'Sulfate'; sw_clrmp_all(ii) = 3;
convert_unit(ii) = 1000*96.06; unit_ti{ii} = 'mg/L';
% ii = ii + 1; plot_comp{ii} = 'S(-2)'; %xlimv(ii,:) = [0 1.5e-4];
% convert_unit(ii) = 1000*32.065; unit_ti{ii} = 'mg/L'; 
% ii = ii + 1; plot_comp{ii} = 'pH'; %xlimv(ii,:) = [-14 14];
ii = ii + 1; plot_comp{ii} = 'FeS(ppt)'; %xlimv(ii,:) = [0 .055]; 
plot_ti{ii} = 'FeS precipitate'; sw_clrmp_all(ii) = 10;
convert_unit(ii) = 1000*(55.845+32.065); unit_ti{ii} = 'mg/Lv';
fl_diff_t1(ii) = 1;
% ii = ii + 1; plot_comp{ii} = 'Siderite'; %xlimv(ii,:) = [0 .055]; 
% fl_diff_t1(ii) = 1;
% ii = ii + 1; plot_comp{ii} = 'Fe(2)'; %xlimv(ii,:) = [0 1];
% % ii = ii + 1; plot_comp{ii} = 'pe'; xlimv(ii,:) = [-14 14];
% % ii = ii + 1; plot_comp{ii} = 'Mn(2)'; xlimv(ii,:) = [-14 14];
% ii = ii + 1; plot_comp{ii} = 'C(4)'; %xlimv(ii,:) = [-14 14];
% ii = ii + 1; plot_comp{ii} = 'C(-4)'; xlimv(ii,:) = [0 0.001];
% ii = ii + 1; plot_comp{ii} = 'O(0)'; fl_diff_t1(ii) = 1;
% ii = ii + 1; plot_comp{ii} = 'Orgcsed'; %fl_diff_t1(ii) = 1;
% ii = ii + 1; plot_comp{ii} = 'Orgcsource'; %fl_diff_t1(ii) = 1;

n_pcomp = ii; plot_comp = plot_comp(1:n_pcomp);

%% ------------------------------------------------------------------------
% load data:
run(fullfile(sim_dir, 'pht3d.m'));

% load time info
times_d = load([sim_dir, slashstr, 'PHT3D_OUTPUT_TIMES.ACN']);

% get domain discr info
fil = fullfile([sim_dir, slashstr, 'pht3d.out']);
fid = fopen(fil, 'r');
while(1) 
    a = fgets(fid);
    if strncmp(a, ' THE TRANSPORT MODEL CONSISTS OF', 32)
        break
    end
end
d = textscan(a, '%s%s%s%s%s%d%s%d%s%d%s'); 
nlay = d{6}; nrow = d{8}; ncol = d{10};
while(1) 
    a = fgets(fid);
    if strncmp(a, '                                       WIDTH ALONG ROWS (DELR)', 62)
        break
    end
end
d = textscan(a, '%s%s%s%s%s%f'); 
DELR = d{6}; % delta along rows [m]
d = textscan(fid, '%s%s%s%s%s%f',1); 
DELC = d{6}; % delta along cols [m]
d = textscan(fid, '%s%s%s%s%s%s%f',1); 
TOP = d{7}; % top elevation [m]
d = textscan(fid, '%s%s%s%s%f%s%s%s', nlay); 
dz_v = d{5}; % vector of layer widths [m]
BOTM = TOP-cumsum(dz_v);
fclose(fid);
z = ([TOP; BOTM(1:end-1)]+BOTM)/2;


% convert to variable name in pht3d.m (check it's there)
plot_comp0 = plot_comp;
for ii = 1: n_pcomp
    plot_comp0{ii}((plot_comp0{ii} == '(')) = '_';
    plot_comp0{ii}(plot_comp0{ii} == '-') = '_';
    plot_comp0{ii}((plot_comp0{ii} == ')')) = ' ';
    plot_comp0{ii} = plot_comp0{ii}(~isspace(plot_comp0{ii}));
    if ~exist(plot_comp0{ii}, 'var')
        fprintf('Error! %s does not exist! Exiting...\n', plot_comp{ii})
        return
    end 
    
end    

if nrow ~= 1, 
    fprintf('Script reads in data for any domain, but only meant for plotting 2D x-sections!! Exiting... \n');
    return
end

% - domain mask and coordinates
% mask = zeros(nlay,ncol);
fid = fopen(ibound_fil);
mask = textscan(fid,'%d',nlay*ncol);
mask = reshape(mask{1}, ncol, nlay)';
xcoord = [DELC: DELC: DELC*double(ncol)];
ycoord_neg = -BOTM;
% return
% plot
fig_i = 1; nfigs = 1;
figure(fig_i)
pagesize = [9 6.5]; % inch, WxH
if fl_mov && fl_recmov
     % separate movie for each figure
     mov = cell(nfigs,1);
     for ff = 1: nfigs
         movfil = ['mov', num2str(ff), '.avi'];
         mov{ff} = VideoWriter(movfil);
         mov{ff}.FrameRate = fps; % fps
         open(mov{ff});
     end
end
      
prow = 3; pcol = 1;
ind = find(sim_dir == '/', 1, 'last');
simname = sim_dir(ind+1:end);
ntimes = length(times_d);
for ts = 1: ntimes
    for ii = 1:n_pcomp
        f_i = ceil(ii/(prow*pcol));
        if f_i > nfigs
            fprintf('Error!  Exceeding nfigs!  Exiting...\n');
            return
        end
        p_i = mod(ii-1,(prow*pcol))+1;
        
        figure(f_i), %orient landscape
        
        % - set window size for movie
        if fl_mov && p_i == 1, 
            set(fig_i, 'PaperUnits', 'inches', ...
                'PaperSize', pagesize);     
            set(fig_i, 'Units', 'inches', 'Position', [0 0 pagesize]); 
            set(gcf, 'Color', [1 1 1]); 
        else
            set(fig_i, 'PaperUnits', 'inches', ...
                'PaperSize', pagesize);
        end  
        
        % - set colorbar
        N_cm0 = 32; % default is 64 bins for colorbar
        switch sw_clrmp_all(ii)
            case 1
                % red with mid yellow
%                 cmcurr0 = cool(N_cm0); %cmcurr0 = cmcurr0(end:-1:1,:); % first range is white
                cmcurr0 = hot(N_cm0); cmcurr0 = cmcurr0(end:-1:1,:); % first range is white
            case 2
                % white to red (thru pink):
                cmcurr0 = nan(N_cm0,3);
                del = 1/(N_cm0-1);
                cmcurr0(:,1) = 1; cmcurr0(:,2:3) = repmat([0:del:1]', 1, 2);
                cmcurr0 = cmcurr0(end:-1:1,:);
            case 3
                % greenscale (some yellow)
                cmcurr0 = hot(N_cm0); cmcurr0 = cmcurr0(end:-1:1,:); % first range is white
                cmcurr0 = cmcurr0(:,[2 1 3]);
            case 4
                % bluescale
                cmcurr0 = hot(N_cm0); cmcurr0 = cmcurr0(end:-1:1,:); % first range is white
                cmcurr0 = cmcurr0(:,[3 2 1]);
            case 5
                % bluegreen
                cmcurr0 = winter(N_cm0-1); cmcurr0 = [1 1 1; cmcurr0(end:-1:1,:)]; % first range is white
            case 6
                % copper
                cmcurr0 = copper(N_cm0-1); cmcurr0 = [1 1 1; cmcurr0(end:-1:1,:)]; % first range is white
            case 7
                % grayscale
                cmcurr0 = gray(N_cm0); 
%                         cmcurr0 = 0.8*cmcurr0(end:-1:1,:); 
                cmcurr0 = 1*cmcurr0(end:-1:1,:); 
                cmcurr0(1,:) = [1 1 1]; % first range is white
            case 8
                % white to green 
                cmcurr0 = nan(n_range,3);
                del = 1/(N_cm0-1);
                cmcurr0(:,2) = 1; cmcurr0(:,[1,3]) = repmat([0:del:1]', 1, 2);
                cmcurr0 = cmcurr0(end:-1:1,:);
            case 0
                % bluegreen
                cmcurr0 = summer(N_cm0-1); cmcurr0 = [1 1 1; cmcurr0(end:-1:1,:)]; % first range is white
            case 10
                % jet
                cmcurr0 = flipud(copper(N_cm0));  
        end        
        
        x = eval(plot_comp0{ii});
        x = reshape(x, ncol, nlay, ntimes);  % ** May need to change order of nlay and ncol!!!
%         ti = [simname, ': ', plot_comp{ii}];
        ti = [plot_ti{ii}];
        x = x * convert_unit(ii);
        ti = [ti, ' ', unit_ti{ii}];
        if fl_diff_t1(ii)
            x = x - repmat(x(:,:,1),[1,1,ntimes]);
            fprintf('Note: %s is difference from t1\n', plot_ti{ii});
%             ti = [ti, ', diff from t1'];
        end        

        if isnan(xlimv(ii,1))
            minval = min(x(:)); 
%             minval = 0;
            xtemp = x(x<1e10); maxval = max(xtemp(:));
            if strcmp(plot_comp{ii}, 'pH'), 
                xtemp = x(x>0);
                minval = min(xtemp(:)); 
            end            
            cax = ([minval, maxval]);
        else
            cax = (xlimv(ii,:));
        end
        minval0 = cax(1) - (cax(2)-cax(1))/N_cm0;
        cax(1) = minval0;
        cmcurr = [1 1 1; cmcurr0];

        % - get data to plot
        xt = x(:,:,ts)'; % data at timestep ts, nlayxncol
        % apply mask
        xt(mask==0) = cax(1);
        
        sp = subplot(prow,pcol,p_i);
        xt(xt>1e10) = nan;
        imagesc(xcoord, ycoord_neg, xt); %colorbar
        hold on
%         if times_d(ts) < 365
%             time_str = [num2str(times_d(ts)), ' days'];
%         else
            yr_str = sprintf('%5.1f', times_d(ts)/365);
            time_str = strtrim([yr_str, ' yrs']);
%         end
        title([ti, ', ', time_str]);
        colormap(sp, cmcurr);
        ylabel('depth [m]'); 
        if p_i == 3
            xlabel('x [m]');
        else
            xlabel('');
        end

        if fl_cbar
            ch = colorbar(sp); 
        end
        %caxis(cax);
        caxis([0 300]);
    end
    drawnow
    
    if fl_mov 
        if fl_recmov
            pause(.1)  % wait n sec 
            MM = getframe(gcf);
            writeVideo(mov{fig_i},MM);
        else
            pause  % wait for press button
        end
    end
end

if fl_mov && fl_recmov, 
    for ff = 1: nfigs
        close(mov{ff}); 
    end
end
if fl_print
    for ff = 1:f_i
        subplot(prow,pcol,2), xlabel('x [m]');
        title('');
        subplot(prow,pcol,1), xlabel('x [m]');
        title('');
        fstr = num2str(ff);
        print('-dtiff', ['-f', fstr], ['fig', fstr, '.tiff']);
    end
end
        
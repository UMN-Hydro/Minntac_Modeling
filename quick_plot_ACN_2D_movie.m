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
sim_dir = '/home/gcng/workspace/ModelRuns_scratch/PHT3D_projects/Minntac/test1/160401a'; % no slash at end
% ************* (end to CUSTOMIZE TO YOUR COMPUTER!!) *********************

fl_print = 1; % 1: print figure to .tiff

% -- what to plot (check pht3d.m in simdir for your options)
n_pcomp_max = 20;
plot_comp = cell(n_pcomp_max,1); 
xlimv = nan(n_pcomp_max, 2);
fl_diff_t1 = zeros(n_pcomp_max,1); 
ii = 0;
ii = ii + 1; plot_comp{ii} = 'Cl'; %xlimv(ii,:) = [0 5];
ii = ii + 1; plot_comp{ii} = 'S(6)'; %xlimv(ii,:) = [0 5];
ii = ii + 1; plot_comp{ii} = 'S(-2)'; %xlimv(ii,:) = [0 1.5e-4];
ii = ii + 1; plot_comp{ii} = 'pH'; %xlimv(ii,:) = [-14 14];
ii = ii + 1; plot_comp{ii} = 'FeS(ppt)'; %xlimv(ii,:) = [0 .055]; 
fl_diff_t1(ii) = 1;
ii = ii + 1; plot_comp{ii} = 'Siderite'; %xlimv(ii,:) = [0 .055]; 
% fl_diff_t1(ii) = 1;
% ii = ii + 1; plot_comp{ii} = 'Fe(2)'; %xlimv(ii,:) = [0 1];
% % ii = ii + 1; plot_comp{ii} = 'pe'; xlimv(ii,:) = [-14 14];
% % ii = ii + 1; plot_comp{ii} = 'Mn(2)'; xlimv(ii,:) = [-14 14];
% ii = ii + 1; plot_comp{ii} = 'C(4)'; %xlimv(ii,:) = [-14 14];
% ii = ii + 1; plot_comp{ii} = 'C(-4)'; xlimv(ii,:) = [0 0.001];
% ii = ii + 1; plot_comp{ii} = 'O(0)'; fl_diff_t1(ii) = 1;
ii = ii + 1; plot_comp{ii} = 'Orgcsed'; %fl_diff_t1(ii) = 1;
ii = ii + 1; plot_comp{ii} = 'Orgcsource'; %fl_diff_t1(ii) = 1;

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


% plot
prow = 2; pcol = 1;
ind = find(sim_dir == '/', 1, 'last');
simname = sim_dir(ind+1:end);
ntimes = length(times_d);
for ts = 1: ntimes
    for ii = 1:n_pcomp
        f_i = ceil(ii/(prow*pcol));
        p_i = mod(ii-1,(prow*pcol))+1;
        figure(f_i), orient landscape
        x = eval(plot_comp0{ii});
        x = reshape(x, ncol, nlay, ntimes);  % ** May need to change order of nlay and ncol!!!
        ti = [simname, ': ', plot_comp{ii}];
        if ~(strcmp(plot_comp0{ii}, 'pH') || strcmp(plot_comp0{ii}, 'pe'))
            x = 1000*x; % mM
            ti = [ti, ' mM'];
        end
        if fl_diff_t1(ii)
            x = x - repmat(x(:,:,1),[1,1,ntimes]);
            ti = [ti, ', diff from t1'];
        end        

        subplot(prow,pcol,p_i)
        xt = x(:,:,ts)'; % data at timestep ts
    %     xt = xt(1:end-1,:);
        xt(xt>1e10) = nan;
        imagesc(xt); colorbar

        if isnan(xlimv(ii,1))
            minval = min(x(:)); 
%             minval = 0;
            xtemp = x(x<1e10); maxval = max(xtemp(:));
            if strcmp(plot_comp{ii}, 'pH'), 
                xtemp = x(x>0);
                minval = min(xtemp(:)); 
            end            
            caxis([minval, maxval]);
        else
            caxis(xlimv(ii,:));
        end
        if times_d(ts) < 365
            time_str = [num2str(times_d(ts)), ' days'];
        else
            yr_str = sprintf('%5.1f', times_d(ts)/365);
            time_str = strtrim([yr_str, ' yrs']);
        end
        title([ti, ', ', time_str]);
    end
    drawnow
    pause
end

% if fl_print
%     for ff = 1:f_i
%         fstr = num2str(ff);
%         print('-dtiff', ['-f', fstr], ['fig', fstr, '.tiff']);
%     end
% end
        

%% fitting both Spatial View and Facing Location
analysis_root = 'C:\Users\dm2574\Dropbox\FMM\data';
monkey = 'Kraut';
% electrode mapping to deuteron channel no.
% load('D:\Dropbox\FMM\data\Kraut\syn_mats\kraut_fmm_grouping_mapping.mat');
% idx = deuteronToGMdriveMapping(:,2)+1;

prs.varname = {'Position',...
               'SpatialView',...
               'SpatialHeading',...
               'EgoBoundary',...
               'Tilt',...          % use the first 2 elements
               'OriSpeed',...      % use the first 2 elements: yaw and pitch
               'Elevation',...
               'TransSpeed',...
               'HeadDir',...
               'LFPhase'
               };

prs.vartype = {'2D', '2D', '2D', '2D', '2D', '2D', '1D', '1D', '1Dcirc', '1Dcirc'};

prs.nbins = {[20,20], [20,20], [20,20], [20,20], [18,6], [16,16], 12, 15, 18, 18};

prs.binrange = {[-1650,-1650;1650,1650], [-5200,-4400;5200,4400], [-5200,-4400;5200,4400],...
                [-1650,-1650;1650,1650], [-0.8,-0.3;1.0,0.3], [-120,-120;120,120],...
                [100,700], [0,1200], [-180,180], [-180,180]};

prs.lambda = {8, 4, 4, 8, 15, 10, 50, 50, 50, 50};     

prs.nfolds = 5;
prs.dt     = 0.02;
prs.filtwidth = 3;
prs.linkfunc  = 'log';
prs.alpha     = 0.05;

fmm_kraut.models_gam_spatialview_spatialheading = cell(size(fmm_kraut,1),1);

for i = 1:size(fmm_kraut,1)
    fprintf(['i = ' num2str(i) '\n']);
    m_d = [monkey, '_', num2str(fmm_kraut.expt_date(i))];
    file_path = fullfile(analysis_root, monkey, m_d, num2str(fmm_kraut.expt_session(i)));
    deuteron = load(fullfile(file_path,'deuteron.mat')); deuteron = deuteron.deuteron;
    behavior = load(fullfile(file_path,'behavior.mat')); behavior = behavior.behavior;
    lfp_path = fullfile('E:\deuteron', m_d, deuteron.expt.session, [m_d '_' deuteron.expt.session '.lfp']);

    ch_lfp = fmm_kraut.chan_deuteron(i);
    lfp = LoadBinary(lfp_path, 'frequency', 1000, 'nChannels', 64, 'channels', ch_lfp);     
    lfp_ph = fmmLFPhase(deuteron,behavior,double(lfp),'passband',[1 10]);

    ts_edges = [behavior.ts, behavior.ts(end)+behavior.ts(2)];
    yt = histcounts(fmm_kraut.spike_ts{i}, ts_edges);
    yt = yt(:);

    nan_ok = isnan(behavior.eye_arena_inter(1,:))|isnan(behavior.head_arena_inter(1,:));
    behavior.position_xy(:,nan_ok) = nan;
    behavior.eye_arena_inter(:,nan_ok) = nan;
    behavior.head_arena_inter(:,nan_ok) = nan;
    behavior.ego_boundary_xy(:,nan_ok) = nan;
    behavior.tilt(:,nan_ok) = nan;
    behavior.yawpitchroll_vel(:,nan_ok) = nan;
    behavior.elevation_z(:,nan_ok) = nan;
    behavior.trans_speed_sm(:,nan_ok) = nan;
    behavior.headdir_azimuth_tilt(:,nan_ok) = nan;
    lfp_ph(:,nan_ok) = nan;

    if isfield(behavior, 'eyepos_xy') > 0.1
        xt = {behavior.position_xy',...
              behavior.eye_arena_inter',...
              behavior.head_arena_inter',...  
              behavior.ego_boundary_xy',...
              behavior.tilt(1:2,:)',...
              medfilt1(behavior.yawpitchroll_vel(1:2,:)',3,[],1),...
              behavior.elevation_z',...
              behavior.trans_speed_sm',...
              behavior.headdir_azimuth_tilt',...
              lfp_ph'};         

        model_this = BuildGAM_Economic_JL_DM(xt, yt, prs);
        fmm_kraut.models_gam_spatialview_spatialheading{i} = model_this;
    end

end

%% remove ith model not fitted (all nans) to save space

for i = 1:size(fmm_kraut,1)
    model_this = fmm_kraut.models_gam_spatialview_spatialheading{i,1};
    if ~isempty(model_this)
        nModel = length(model_this.class);
        nan_model = [];
        for j = 1:nModel
            if isnan(nanstd(model_this.testFit{j,1}(:)))
                nan_model = cat(1,nan_model,j);
            end
        end
        
        model_this.testFit(nan_model)  = {[]};
        model_this.trainFit(nan_model) = {[]};
        model_this.wts(nan_model)      = {[]};
        model_this.marginaltunings(nan_model) = {[]};
        
        fmm_kraut.models_gam_spatialview_spatialheading{i,1} = model_this;
    end
end

%% fitting both Head Dir. and Gaze Dir.
analysis_root = 'C:\Users\dm2574\Dropbox\FMM\data';
monkey = 'Kraut';

prs.varname = {'Position',...
               'SpatialHeading',...
               'EgoBoundary',...
               'Tilt',...          % use the first 2 elements
               'OriSpeed',...      % use the first 2 elements: yaw and pitch
               'Elevation',...
               'TransSpeed',...
               'HeadDir',...
               'GazeDir',...
               'LFPhase'
               };

prs.vartype = {'2D', '2D', '2D', '2D', '2D', '1D', '1D', '1Dcirc', '1Dcirc', '1Dcirc'};

prs.nbins = {[20,20], [20,20], [20,20], [18,6], [16,16], 12, 15, 18, 18, 18};

prs.binrange = {[-1650,-1650;1650,1650], [-5200,-4400;5200,4400],...
                [-1650,-1650;1650,1650], [-0.8,-0.3;1.0,0.3], [-120,-120;120,120],...
                [100,700], [0,1200], [-180,180], [-180,180], [-180,180]};

prs.lambda = {8, 4, 8, 15, 10, 50, 50, 50, 50, 50};     

prs.nfolds = 5;
prs.dt     = 0.02;
prs.filtwidth = 3;
prs.linkfunc  = 'log';
prs.alpha     = 0.05;

fmm_kraut.models_gam_headdir_eyedir = cell(size(fmm_kraut,1),1);

for i = 1:size(fmm_kraut,1)
    fprintf(['i = ' num2str(i) '\n']);

    m_d = [monkey, '_', num2str(fmm_kraut.expt_date(i))];
    file_path = fullfile(analysis_root, monkey, m_d, num2str(fmm_kraut.expt_session(i)));
    deuteron = load(fullfile(file_path,'deuteron.mat')); deuteron = deuteron.deuteron;
    behavior = load(fullfile(file_path,'behavior.mat')); behavior = behavior.behavior;
    lfp_path = fullfile('D:\FMM\deuteron', m_d, deuteron.expt.session, [m_d '_' deuteron.expt.session '.lfp']);

    ch_lfp = fmm_kraut.chan_deuteron(i);
    lfp = LoadBinary(lfp_path, 'frequency', 1000, 'nChannels', 64, 'channels', ch_lfp);     
    lfp_ph = fmmLFPhase(deuteron,behavior,double(lfp),'passband',[1 10]);

    ts_edges = [behavior.ts, behavior.ts(end)+behavior.ts(2)];
    yt = histcounts(fmm_kraut.spike_ts{i}, ts_edges);
    yt = yt(:);

    nan_ok = isnan(behavior.eye_arena_inter(1,:))|isnan(behavior.head_arena_inter(1,:));
    behavior.position_xy(:,nan_ok) = nan;
    behavior.head_arena_inter(:,nan_ok) = nan;
    behavior.ego_boundary_xy(:,nan_ok) = nan;
    behavior.tilt(:,nan_ok) = nan;
    behavior.yawpitchroll_vel(:,nan_ok) = nan;
    behavior.elevation_z(:,nan_ok) = nan;
    behavior.trans_speed_sm(:,nan_ok) = nan;
    behavior.headdir_azimuth_tilt(:,nan_ok) = nan;
    behavior.eyedir_azimuth_earth(:,nan_ok) = nan;
    lfp_ph(:,nan_ok) = nan;

    if isfield(behavior, 'eyepos_xy') > 0.1
        xt = {behavior.position_xy',...
              behavior.head_arena_inter',...  
              behavior.ego_boundary_xy',...
              behavior.tilt(1:2,:)',...
              medfilt1(behavior.yawpitchroll_vel(1:2,:)',3,[],1),...
              behavior.elevation_z',...
              behavior.trans_speed_sm',...
              behavior.headdir_azimuth_tilt',...
              behavior.eyedir_azimuth_earth',...
              lfp_ph'};         

        model_this = BuildGAM_Economic_JL_DM(xt, yt, prs);
        fmm_kraut.models_gam_headdir_eyedir{i} = model_this;
    end
 
end

%% remove ith model not fitted (all nans) to save space

for i = 1:size(fmm_kraut,1)
    model_this = fmm_kraut.models_gam_headdir_eyedir{i,1};
    if ~isempty(model_this)
        nModel = length(model_this.class);
        nan_model = [];
        for j = 1:nModel
            if isnan(nanstd(model_this.testFit{j,1}(:)))
                nan_model = cat(1,nan_model,j);
            end
        end
        
        model_this.testFit(nan_model)  = {[]};
        model_this.trainFit(nan_model) = {[]};
        model_this.wts(nan_model)      = {[]};
        model_this.marginaltunings(nan_model) = {[]};
        
        fmm_kraut.models_gam_headdir_eyedir{i,1} = model_this;
    end
end

%% fitting just eye position and eye velocity in head
analysis_root = 'D:\Dropbox\FMM\data';
monkey = 'Kraut';

prs.varname = {'EyePos',...
               'EyeVel'
               };

prs.vartype = {'2D', '2D'};

prs.nbins = {[20,20], [20,20]};

prs.binrange = {[-40,-30;40,30], [-100,-100;100,100]};

prs.lambda = {10, 10};

prs.nfolds = 5;
prs.dt     = 0.02;
prs.filtwidth = 3;
prs.linkfunc  = 'log';
prs.alpha     = 0.05;

fmm_kraut.models_eyeinhead = cell(size(fmm_kraut,1),1);

for i = 1:size(fmm_kraut,1)
    fprintf(['i = ' num2str(i) '\n']);

    m_d = [monkey, '_', num2str(fmm_kraut.expt_date(i))];
    file_path = fullfile(analysis_root, monkey, m_d, num2str(fmm_kraut.expt_session(i)));
    deuteron = load(fullfile(file_path,'deuteron.mat'));  deuteron = deuteron.deuteron;
    behavior = load(fullfile(file_path,'behavior.mat'));  behavior = behavior.behavior;
    ts_edges = [behavior.ts, behavior.ts(end)+behavior.ts(2)];
    yt = histcounts(fmm_kraut.spike_ts{i}, ts_edges);
    yt = yt(:);

    if isfield(behavior, 'eyepos_xy') > 0.1
        eye_pos = behavior.eyepos_xy';
%         eye_pos = fillmissing(eye_pos, 'movmedian', 8, 1, 'EndValues', 'none');
        eye_pos_ds = JL_Desaccader_Simple(eye_pos);
        eye_vel = diff(eye_pos_ds,1,1) ./ 0.02;  % deg/s
        eye_vel = [eye_vel(1,:); eye_vel];
        eye_vel = smoothdata(eye_vel, 1, 'movmedian', 8, 'omitnan');
        xt = {eye_pos,...
              eye_vel
              };         

        model_this = BuildGAM_Economic_JL_DM(xt, yt, prs);
        fmm_kraut.models_eyeinhead{i} = model_this;
    end

end

%%
for i = 1:size(fmm_kraut,1)
    model_this = fmm_kraut.models_eyeinhead{i,1};
    if ~isempty(model_this)
        nModel = length(model_this.class);
        nan_model = [];
        for j = 1:nModel
            if isnan(nanstd(model_this.testFit{j,1}(:)))
                nan_model = cat(1,nan_model,j);
            end
        end
        
        model_this.testFit(nan_model)  = {[]};
        model_this.trainFit(nan_model) = {[]};
        model_this.wts(nan_model)      = {[]};
        model_this.marginaltunings(nan_model) = {[]};
        
        fmm_kraut.models_eyeinhead{i,1} = model_this;
    end
end

%%
load('C:\Users\dm2574\Dropbox\FMM\data\Kraut\syn_mats\archive\kraut_fmm_headVSeye_eyeinhead.mat');

nVar    = 2;
nNeuron = size(fmm_kraut,1);
sig_var = zeros(nNeuron,nVar);
codedvar_bm = zeros(1,nVar);
testfit_mn  = nan(nNeuron,nVar);

for i = 1:nNeuron
    model = fmm_kraut.models_eyeinhead{i};
    if ~isempty(model)
        for iVar = 1:nVar
            tf_this = model.testFit{iVar}(:,3);
            p = signrank(tf_this, 0, 'Tail', 'right');
            if p < 0.05
                sig_var(i,iVar) = 1;
                testfit_mn(i,iVar) = mean(tf_this);
            end
        end
        bm = model.bestmodel;
        tf_this = model.testFit{bm}(:,3);
        p = signrank(tf_this, 0, 'Tail', 'right');
        if p < 0.05
            codedvar_bm = codedvar_bm + model.class{bm};
        end
    end
end

%% example eye pos and eye vel. tuning curve
idx_pos = find(sig_var(:,1));
idx_vel = find(sig_var(:,2));
i = 21; i = idx_pos(i);
figure
subplot(1,2,1)
zz = fmm_kraut.models_eyeinhead{i,1}.marginaltunings{1,1}{1,1};
tt = reshape(fmm_kraut.models_eyeinhead{i,1}.nDataPoints{1,1},prs_eyeinhead.nbins{1})*0.02;
c_lim = round([0 max(zz(:))],2);
b = imagesc(zz,c_lim); box off; colorbar; pbaspect([4 3 1]);
set(b, 'AlphaData', ~(tt<0.3));
% set(gca,'YDir', 'normal');
set(gca, 'XTick', 0.5:5:20.5, 'XTickLabel',-40:20:40);
set(gca, 'YTick', 0.5:5:20.5, 'YTickLabel',30:-15:-30);
xlabel('Hor. (°)'); ylabel('Ver. (°)');

i = 9; i = idx_vel(i);
subplot(1,2,2)
zz = fmm_kraut.models_eyeinhead{i,1}.marginaltunings{1,2}{1,2};
tt = reshape(fmm_kraut.models_eyeinhead{i,1}.nDataPoints{1,2},prs_eyeinhead.nbins{2})*0.02;
c_lim = round([0 max(zz(:))],2);
b = imagesc(zz,c_lim); axis square; colorbar; box off;
set(b, 'AlphaData', ~(tt<0.3));
set(gca, 'XTick', 0.5:5:20.5, 'XTickLabel', -150:75:150);
set(gca, 'YTick', 0.5:5:20.5, 'YTickLabel', 150:-75:-150);
xlabel('Hor. Vel. (°/s)'); ylabel('Ver. Vel. (°/s)');

%% significant saccade-related responses
analysis_root = 'D:\Dropbox\FMM\data';
sac_sig = zeros(size(fmm_kraut,1),3);
for i = 1:size(fmm_kraut,1)
    monkey = fmm_kraut.animal_id{i};
    m_d    = [monkey, '_', num2str(fmm_kraut.expt_date(i))];
    file_path = fullfile(analysis_root, monkey, m_d, num2str(fmm_kraut.expt_session(i)));
    behavior  = load(fullfile(file_path,'behavior.mat')); behavior = behavior.behavior;
    ts_spk    = fmm_kraut.spike_ts{i};
    ts = behavior.ts';
    eye_pos = behavior.eyepos_xy';
    [~, saccade] = JL_Desaccader_Simple(eye_pos);
    ts_sac = ts(saccade.sac_st);
    [~, idx_pre] = Sync(ts_spk, ts_sac, 'durations', [-0.4 0]);
    [~,  idx_in] = Sync(ts_spk, ts_sac, 'durations', [0 0.4]);
    ct_pre = []; ct_in = [];
    for j = 1:length(ts_sac)
        ct_pre = cat(1, ct_pre, sum(idx_pre==j));
        ct_in  = cat(1, ct_in,  sum(idx_in==j));
    end
    if ttest(ct_pre, ct_in, 'Alpha', 0.01)
        sac_sig(i,1) = 1;
        if mean(ct_pre) > mean(ct_in)
            sac_sig(i,2) = 1;
        else
            sac_sig(i,3) = 1;
        end
    end
end

%% saccade aligned responses
idx_sig = find(sac_sig(:,1));
i = 41;
i = idx_sig(i);
monkey = fmm_kraut.animal_id{i};
m_d = [monkey, '_', num2str(fmm_kraut.expt_date(i))];
file_path = fullfile(analysis_root, monkey, m_d, num2str(fmm_kraut.expt_session(i)));
behavior = load(fullfile(file_path,'behavior.mat')); behavior = behavior.behavior;
ts_spk = fmm_kraut.spike_ts{i};
ts = behavior.ts';
eye_pos = behavior.eyepos_xy';
[eye_pos_ds, saccade] = JL_Desaccader_Simple(eye_pos);
ts_sac = ts(saccade.sac_st);

[synced, idx] = Sync(ts_spk, ts_sac, 'durations', [-0.4 0.8]);
[mn, err, ts_bins] = SyncHist(synced, idx, 'mode', 'mean', 'durations', [-0.4 0.8], 'nBins', 20, 'error', 'sem');
errorbar(ts_bins, mn, err); box off
xlabel('Time relative to saccade onset (s)');
ylabel('Firing rate (Hz)');
xlim([-0.45 0.85]);

%% example eye movement plot, shaded saccade period
figure
ok = ts > 16 & ts < 26;
plot(ts(ok), eye_pos(ok,1),'LineWidth',2); hold on
plot(ts(ok), eye_pos(ok,2),'LineWidth',2); box off
% legend('Hor.','Ver.');
ylim([-40 40]); ylabel('°');
xx = ts(saccade.sac_st);
ok = xx > 16 & xx < 26;
xx = xx(ok);
plot([xx xx], [-40 40], '--', 'Color', [0 0 0]+0.5);

%% bar plot 
ok_hpc = strcmp(fmm_kraut.region,'HPC') & ~isempty(fmm_kraut.models_eyeinhead);
eyepos_hpc = sum(sig_var(ok_hpc,1))/sum(ok_hpc);
eyevel_hpc = sum(sig_var(ok_hpc,2))/sum(ok_hpc);
% sac_hpc = sum(sac_sig(ok_hpc,1))/sum(ok_hpc);

ok_ec = strcmp(fmm_kraut.region,'EC') & ~isempty(fmm_kraut.models_eyeinhead);
eyepos_ec = sum(sig_var(ok_ec,1))/sum(ok_ec);
eyevel_ec = sum(sig_var(ok_ec,2))/sum(ok_ec);
% sac_ec = sum(sac_sig(ok_ec,1))/sum(ok_ec);

ok_sub = ~(strcmp(fmm_kraut.region,'EC')|strcmp(fmm_kraut.region,'HPC')) & ~isempty(fmm_kraut.models_eyeinhead);
eyepos_sub = sum(sig_var(ok_sub,1))/sum(ok_sub);
eyevel_sub = sum(sig_var(ok_sub,2))/sum(ok_sub);
% sac_sub = sum(sac_sig(ok_sub,1))/sum(ok_sub);

y = [eyepos_hpc, eyepos_ec, eyepos_sub;...
     eyevel_hpc, eyevel_ec, eyevel_sub;...
%      sac_hpc, sac_ec, sac_sub];
    ];

figure
bar(y); box off
set(gca, 'XTick', 1:3, 'XTickLabel', {'Eye Pos.', 'Eye Vel.', 'Saccade'});
xlabel('Encoding variable'); ylabel('Fraction');
legend({['HPC (n=' num2str(sum(ok_hpc)) ')'],...
        ['EC (n=' num2str(sum(ok_ec)) ')'],...
        ['SUB (n=' num2str(sum(ok_sub)) ')']});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% tuning curve based approach
i = 36;
fprintf(['i = ' num2str(i) '\n']);
m_d = [monkey, '_', num2str(fmm_kraut.expt_date(i))];
file_path = fullfile(analysis_root, monkey, m_d, num2str(fmm_kraut.expt_session(i)));
deuteron = load(fullfile(file_path,'deuteron.mat')); deuteron = deuteron.deuteron;
behavior = load(fullfile(file_path,'behavior.mat')); behavior = behavior.behavior;
ts_spk = fmm_kraut.spike_ts{i};
ts = behavior.ts';

eye_pos = behavior.eyepos_xy';
eye_pos = JL_Desaccader_Simple(eye_pos);
eye_vel = diff(eye_pos,1,1);
eye_vel = [eye_vel(1,:); eye_vel];
eye_vel = smoothdata(eye_vel, 1, 'movmedian', 8, 'includenan');

eye_pos(:,1) = (eye_pos(:,1)+30)/60;
eye_pos(:,2) = (eye_pos(:,2)+25)/50;
eye_pos(eye_pos<0 | eye_pos>1) = nan;
fm = Map([ts, eye_pos], ts_spk, 'nBins', [30 30], 'Smooth', [0 0]);
fm.z(fm.time<0.5) = nan;
fm.time(fm.time<0.5) = nan;
c_lim  = round([nanmin(fm.z(:)) nanmax(fm.z(:))],2);

subplot(1,2,1)
b = imagesc(fm.z, c_lim); axis square; axis square
set(b, 'AlphaData', ~(isnan(fm.time)));
set(gca, 'YDir', 'normal');

eye_vel(:,1) = (eye_vel(:,1)+0.8)/1.6;
eye_vel(:,2) = (eye_vel(:,2)+0.8)/1.6;
eye_vel(eye_vel<0 | eye_vel>1) = nan;
fm = Map([ts, eye_vel], ts_spk, 'nBins', [30 30], 'Smooth', [0 0]);
fm.z(fm.time<0.5) = nan;
fm.time(fm.time<0.5) = nan;
c_lim  = round([nanmin(fm.z(:)) nanmax(fm.z(:))],2);

subplot(1,2,2)
b = imagesc(fm.z, c_lim); axis square; axis square
set(b, 'AlphaData', ~(isnan(fm.time)));
set(gca, 'YDir', 'normal');

%% tuning curve shuffle analysis
[info_eyepos_data, info_eyepos_shfl, info_eyevel_data, info_eyevel_shfl] = fmmEyeinHeadAnalysis(fmm_kraut);

%%
[expt_uq,ia,ic]  = unique(fmm(:,1:3),'rows');
[gc,ga] = groupcounts(ic);

%%
i = 16;
expt = expt_uq(i,:);
ts_spk = fmm.spike_ts(ic==i,:);
elec_clu = fmm(ic==i,5:6); elec_clu = table2array(elec_clu);
region = fmm.region(ic==i,:);
md_spk = fmm.models_gam_chunk(ic==i,:);
fmmdirs  = fmmGetDir(expt);
[deuteron, behavior, tetrode] = fmmLoadData(fmmdirs);
session_flag = fmm.session_flag(ic==i,:);

%
for j = 1:30
    if strcmp(session_flag{j},'good')
        clu = elec_clu(j,:);
        if strcmp(expt.animal_id,'Lysander')&&(expt.expt_date<= 20181030 && expt.expt_date > 20180418)
            clu(:,1) = 9-clu(:,1);
        end
        fmmPlotRawTuningCurve(behavior, ts_spk(j), md_spk(j), tetrode, clu);
        sgtitle([expt.animal_id{1} ' ' num2str(expt.expt_date) ' clu ' num2str(clu)]);
    end
end

%%
llh = fmmTableToLLHQuant(fmm);

%% tuning curve, data and model-fit, 1st order model
i = 287;
fmmPlotTuningCurveDataModel(fmm(i,:));

%% tuning curve, data and model-fit, best model
fmmPlotTuningCurveDataModelBest(fmm(i,:));

%% example tuning curve and LLH plots
egs = [89,238,287,322,384,404,441,476,597];
ii = 10;
i = egs(ii);
fmmPlotTuningCurveDataModel(fmm(i,:));

% log likelihood plots, all models
fmmPlotLLH(fmm(i,:));

%% spatial distribution of preferred 'receptive field', cluster analysis using PCA
fmmPlotDistrbRFcluster(fmm);

%% example SH, Tt, Pos, and EB tuning curves, plot 20 each
fmmPlotManyEgTCs(fmm);

%% additive or multiplicative, compare fraction of explained variance
% fmmPlotEVaddVSmultip(fmm);

%% plot the most co-coded variables across regions
fmmPlotModelTypeVarCo(fmm);

%% farction of encoding neurons across regions, individual vars, classifications
fmmPlotEncodingFraction(fmm);
fmmPlotEncodingFractionHPCsubregions;

%% fraction of encoding neurons vs. sorting quality
fmmPlotEncodingFractionSortingQuality(fmm);

%% fraction of encoding neurons across head-body-tail HPC, MEC vs. other EC
fmmPlotEncodingFractionAnatDstrb(fmm);

%% firing rate distribution across regions
fmmPlotFiringRateLog(fmm);

%% grid cell analysis
fmmGridCellAnalysis(fmm);
fmmPlotGridTuningAnalysis;

%% place cell analysis
fmmPlaceCellAnalysis(fmm);

%% relationship between tuning and putative cell types
fmmTuningWaveform;

%% traditional speed tuning and head direction analysis
fmmSpeedAnalysis(fmm);
fmmHeadDirAnalysis(fmm);

%% summary of tranditional analyses
fmm_traditional_analysis;

%% relationship between coding of different variables
% the effect of gaze on the sharpness of tuning
fmmCorrGeometryVarsGazeUncertainty;

%% Lysander: task vs. no task
fmmLysanderTaskVSNoTask;

%% 1. behavioral var. distribution
[expt_uq,ia,ic]  = unique(fmm(:,1:3),'rows');
nCount = cell(length(ia),1);
for i = 1:length(ia)
    expt_this = expt_uq(i,:);
    fmmdirs  = fmmGetDir(expt_this);
    [~, behavior, ~] = fmmLoadData(fmmdirs);
    [xc, nCount(i,:)] = fmmBehaviorDstrb(behavior);
end

%% 1.1 eye-in-head distribution for Kraut
ok = strcmp(fmm.animal_id,'Kraut');
[expt_uq,ia,ic]  = unique(fmm(ok,1:3),'rows');
nCount = [];
for i = 1:length(ia)
    expt_this = expt_uq(i,:);
    fmmdirs  = fmmGetDir(expt_this);
    [~, behavior, ~] = fmmLoadData(fmmdirs);
    if isfield(behavior,'eyepos_xy')
        [xc, nCount_this] = fmmBehaviorDstrbEyeInHead(behavior);
        nCount = cat(1, nCount, nCount_this);
    end
end

%% 2. behavioral var. distribution
ok = strcmp(expt_uq.animal_id,'Bruno');
fmmPlotBehaviorDstrb(xc, nCount(ok,:));

%% 3. behavioral var. distribution, position bin count distribution
nCount_pos = nCount(:,1); nCount_pos = cell2mat(nCount_pos);
for k = 1:136; for i=1:20; for j =1:20; if sqrt(xc{1}{1,1}(i).^2+xc{1}{1,2}(i).^2) > 1650; idx = sub2ind([20,20],i,j); nCount_pos(k,idx)=1; end; end; end; end
nCount_pos_empty = sum(nCount_pos==0,2);

%% control analysis: spatial view vs. spatial heading; head azimuth dir. vs. eye azimuth dir.
fmmPlotControlHeadVSEye(fmm_kraut,behavior);

%% LFP analysis, power spectrum
i = 200;
fmmLFPSpectralFooof(fmm(i,:));

%% LFP analysis, plot wavelet-based spectrogram
i = 200;
fmmPlotLFPSpectrogramCWT(fmm(i,:));

%% LFP analysis, average spectral density over sessions, individual regions, individual monkeys
fmmLFPSpectrogramCWT(fmm);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% temp
bad_sum = zeros(1,599);
for i = 1:599
    ts_this = fmm.spike_ts{i};
    not_ok  = find(diff(ts_this) < 0.002)+1;
    bad_sum(i) = length(not_ok);
    ts_this(not_ok) = [];
    fmm.spike_ts{i} = ts_this;
end

%%
i = 598;
ts_this = fmm.spike_ts{i};
[ccg,t] = CCG(ts_this,ones(size(ts_this))); 
bar(t,ccg,'k','EdgeColor','k'); axis square; box off;
xlabel('dt (s)'); ylabel('Count');
    
    

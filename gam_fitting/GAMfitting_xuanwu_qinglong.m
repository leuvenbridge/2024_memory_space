clear; clc
monkey = 'qinglong';
rootData   = fullfile('C:\Users\pc\OneDrive\MyLab\data', monkey);

sessions = {{'20230309'},{'20230310'},{'20230313'},{'20230315'},...
            {'20230316'},{'20230320'},{'20230322'},{'20230323'},{'20230324'},...
            {'20230327'},{'20230328'},{'20230329'},{'20230331'},{'20230404'},...
            {'20230406'},{'20230410'},{'20230411'},{'20230412'},{'20230414'},...
            {'20230417'},{'20230418'},{'20230420'},{'20230425'},{'20230428'},...
            {'20230505'},{'20230508'},{'20230512'},{'20230516'},{'20230519'},...
            {'20230523'},{'20230525'},{'20230531'},{'20230602'},{'20230606'},...
            {'20230607'},{'20230609'}};

%%
clear; clc
monkey = 'xuanwu';
rootData = fullfile('C:\Users\pc\OneDrive\MyLab\data', monkey);
sessions = {{'20230216'},{'20230219'},{'20230221'},{'20230222'},{'20230224'},...
            {'20230225'},{'20230228'},{'20230301'},{'20230303'},{'20230305'},...
            {'20230307'},{'20230308'},{'20230309'},{'20230311'},{'20230320'},...
            {'20230321'},{'20230327'},{'20230329'},{'20230331'},{'20230401'},...
            {'20230403'},{'20230405'},{'20230410'},{'20230417'},{'20230419'},...
            {'20230421'},{'20230422'},{'20230424'},{'20230426'},{'20230504'},...
            {'20230506'},{'20230508'},{'20230510'},{'20230531'},{'20230602'},...
            {'20230720'},{'20230725'},{'20230801'},{'20230802'},{'20230804'},...
            {'20230807'},{'20230830'},{'20230914'},{'20230915'},{'20230918'},...
            {'20230921'},{'20230927'}};

%%
% downsample factor, 0.02 s step is fine in general
% prs.dt needs to be changed accordingly if ds_factors is changed
ds_factor = 2;
% names of all variables you want to fit
prs.varname = {'Position',...       % self location
               'SpatialView',...    % spatial view
               'EgoCenter',...      % egocentric center
               'Tilt',...           % head tilt, use the first 2 elements: pitch and roll
               'OriSpeed',...       % head orientation speed, use the first 2 elements: yaw and pitch
               'Elevation',...      % head height
               'TransSpeed',...     % translational speed
               'HeadDir',...        % horizontal head direction
               };

prs.vartype = {'2D', '2D', '2D', '2D', '2D', '1D', '1D', '1Dcirc'};

prs.nbins = {[20,20], [20,20], [20,20], [18,9], [16,16], 15, 15, 15};

prs.binrange = {[-1750,-1750;1750,1750],...
                [-7000,-4600;7000,4600], [-2400,-2400;2400,2400],...
                [-0.6,-0.3;0.6,0.3], [-120,-120;120,120],...
                [100,1000], [0,1500], [-180,180]};

prs.lambda = {8, 4, 8, 15, 10, 50, 50, 50};     

prs.nfolds    = 5;             % n folds cross-validation
prs.dt        = 0.02;
prs.filtwidth = 3;
prs.linkfunc  = 'log';
prs.alpha     = 0.05;

nSession = length(sessions);

for iSession = 1:nSession
    session_this   = sessions{iSession}{1}
    session_date   = strsplit(session_this, '_');
    session_date   = session_date{1};
    datapath       = fullfile(rootData, session_date, 'ff');
    modelpath      = fullfile(datapath, 'models_gam_spykingcircus.mat');
    spikespath     = fullfile(datapath, 'spikes_spykingcircus.mat');
%     if exist(spikespath,'file') && isfolder(datapath) && ~exist(modelpath,'file')
    
    % temporary: 20231202
    datestr = '30-11-2023';
    formatin = 'dd-mm-yyyy';
    olddatenum = datenum(datestr, formatin);
    if exist(spikespath,'file') && dir(spikespath).datenum > olddatenum
    % temporary: 20231202
    
        syncpath       = fullfile(datapath, 'sync.mat');      load(syncpath);
        spikespath     = fullfile(datapath, 'spikes_spykingcircus.mat');    load(spikespath);
        viconpath      = fullfile(datapath, 'vicon.mat');     load(viconpath);
        eyecalib       = fullfile(datapath, 'binocalib.mat'); load(eyecalib);
        ffeye          = fullfile(datapath, 'ffbinoeye.mat'); load(ffeye);
        cd(datapath);
        if strcmp(monkey, 'xuanwu')
            eyecalibpara   = GetBinoEyeCalib_xuanwu(BinoCalib);
            behavior       = ExtractBehavior_xuanwu(vicon, FFEye, eyecalibpara, sync);

        elseif strcmp(monkey, 'qinglong')
            eyecalibpara   = GetBinoEyeCalib_qinglong(BinoCalib);
            behavior       = ExtractBehavior_qinglong(vicon, FFEye, eyecalibpara, sync);
        end

        behavior.yawpitchroll_vel = medfilt1(behavior.yawpitchroll_vel,10,[],2);

        ts = behavior.ts;
        ts = ts(1:ds_factor:end);
        dt = ts(2) - ts(1);
        ts_edges = [ts, ts(end)+dt];

        nUnit = size(spikes,1);
        models_gam = cell(nUnit,1);

        xt = {behavior.position_xy(:,1:ds_factor:end)',...
              behavior.eye_left_arena_inter(:,1:ds_factor:end)',... 
              behavior.ego_center_xy(:,1:ds_factor:end)',...
              behavior.tilt(1:2,1:ds_factor:end)',...
              behavior.yawpitchroll_vel(1:2,1:ds_factor:end)',...
              behavior.elevation_z(:,1:ds_factor:end)',...
              behavior.trans_speed_sm(:,1:ds_factor:end)',...
              behavior.headdir_azimuth_tilt(:,1:ds_factor:end)'};  

        for i = 1:nUnit
            fprintf(['iUnit = ' num2str(i) '\n']);
            yt = histcounts(double(spikes.times{i,1} + sync.SG_first_ts)./sync.SG_samplerate, ts_edges);
            yt = yt(:);
            model_this = BuildGAM_Economic(xt, yt, prs);        

            % remove ith model not fitted (all nans) to save space
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

                models_gam{i,1} = model_this;
            end
        end
        save('models_gam_spykingcircus.mat','models_gam');
    end
end

%% visualize fitting results
% PlotGAM(models_gam{i}, prs);
ffPlotModelTypeVarCo;
ffPlotLLH;

% NOTE: you will need to modify the following 2 scripts to visualize actual
% and fitted tuning curves
ffPlotTuningCurveDataModel;
ffPlotTuningCurveDataModelBest;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
[expt_uq,ia,ic]  = unique(models_gam(:,1:3),'rows');
[gc,ga] = groupcounts(ic);

ic_ok = ga(gc>=9);
expt_uq_ok = expt_uq(ic_ok,:);


%% fitting GAM model with spike history and coupling


%%
clear; clc
monkey = 'xuanwu';
rootData   = fullfile('C:\Users\pc\OneDrive\MyLab\data', monkey);
sessions = {{'20230216'},...
            {'20230221'},...
            {'20230222'},...
            {'20230224'},...
            {'20230225'},...
            {'20230228'},...
            {'20230301'},...
            {'20230303'},...
            {'20230305'},...
            {'20230307'},...
            {'20230308'},...
            {'20230309'},...
            {'20230311'},...
            {'20230320'},...
            {'20230321'},...
            {'20230327'},...
            {'20230401'},...
            {'20230403'},...
            {'20230419'},...
            {'20230417'},...
            {'20230419'},...
            {'20230421'},...
            {'20230422'},...
            {'20230424'},...
            {'20230426'},...
            {'20230504'},...
            {'20230506'},...
            {'20230508'},...
            {'20230510'}};

%%
clear; clc
monkey = 'qinglong';
rootData   = fullfile('C:\Users\pc\OneDrive\MyLab\data', monkey);
sessions = {{'20230404'},...
            {'20230406'},...
            {'20230411'},...
            {'20230412'},...
            {'20230414'},...
            {'20230417'},...
            {'20230418'},...
            {'20230420'},...
            {'20230428'},...
            {'20230505'},...
            {'20230508'},...
            {'20230512'},...
            {'20230516'},...
            {'20230519'},...
            {'20230523'},...
            {'20230525'}};
        
%%
date = {}; group = []; channel = []; cluster = [];
label = {}; fr = []; tp_dur = []; models = {};
spikes_ff   = table(date, group, channel, cluster, label, fr, tp_dur, models);

ii = 0;

for iSession = 2%1:length(sessions)
    
    session_this   = sessions{iSession}{1};
    session_date   = strsplit(session_this, '_');
    session_date   = session_date{1};
    
    datapath       = fullfile(rootData, session_date, 'ff');
    syncpath       = fullfile(datapath, 'sync.mat');      load(syncpath);
    spikespath     = fullfile(datapath, 'spikes.mat');        load(spikespath);
    modelpath      = fullfile(datapath, 'models_gam.mat');    load(modelpath);
    
    nUnit          = size(spikes, 1);
    
    for iUnit = 1:nUnit
        ii = ii + 1;
        ts_spk   = double(spikes.times{iUnit,1} + sync.SG_first_ts)/sync.SG_samplerate;
        wv_spk   = mean(double(spikes.waveform{iUnit,1}),2,'omitnan');
        [~,idx1] = min(wv_spk);
        [~,idx2] = max(wv_spk(idx1:end));
        tp_dur   = idx2 ./ 30;   % trough to peak duration, in ms
        
        spikes_ff.date{ii}      = session_date;
        spikes_ff.group(ii)     = spikes.group{iUnit};
        spikes_ff.channel(ii)   = spikes.channel{iUnit};
        spikes_ff.cluster(ii)   = spikes.cluster{iUnit};
        spikes_ff.label{ii}     = spikes.label{iUnit};
        spikes_ff.fr(ii)        = length(ts_spk)/(ts_spk(end)-ts_spk(1));
        spikes_ff.tp_dur(ii)    = tp_dur;
        spikes_ff.models{ii}    = models_gam{iUnit};
    end

end

%% visualize fitting results
% PlotGAM(models_gam{i}, prs);
ffPlotModelTypeVarCo;
ffPlotLLH;

% NOTE: you will need to modify the following 2 scripts to visualize actual
% and fitted tuning curves
ffPlotTuningCurveDataModel;
ffPlotTuningCurveDataModelBest;












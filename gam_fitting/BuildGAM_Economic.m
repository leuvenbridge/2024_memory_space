function models = BuildGAM_Economic(xt,yt,prs)

% BUILDGAM: Fit generalised additive models to neural data.
%   Model: g(r) = sum_i[f_i(x_i)], where r denotes the mean response, x_i 
%   is the i-th input variable and f_i(x_i) is the tuning of y to x_i
%   (fit to data), and g is the link function denoting neuronal 
%   nonlinearity. Several variants of the model are fit, where each variant
%   corresponds to a unique combination of the input variables x_i. The
%   variant with the highest likelihood is chosen via forward search.
%
% INPUTS:
% xt is an 1 x N cell array where the i-th cell is an array containing 
% values of the i-th variable x_i. Each array is either T x 1 or T x 2 
% depending on whether x_i is 1-dimensional (e.g. speed) or 2-dimensional 
% (e.g. position), where T is the total number of observations.
%
% yt is a T x 1 array of spike counts
%
% prs is a structure specifying analysis parameters. The contents of the
% structure are as follows and must be entered in this order.
% prs.varname   : 1 x N cell array of names of the input variables.
%                 only used for labeling plots
% prs.vartype   : 1 x N cell array of types ('1D','1Dcirc' or '2D') of the input variables. 
%                 used for applying smoothness penalty on tuning functions
% prs.nbins     : 1 x N cell array of number of bins to discretise input variables. 
%                 determines the resolution of the tuning curves
% prs.binrange  : 1 x N cell array of 2 x 1 vectors specifying lower and upper bounds of input variables. 
%                 used to determine bin edges
% prs.nfolds    : Number of folds for cross-validation.
% prs.dt        : Time between consecutive observation samples. 
%                 used for converting weights f_i to firing rate
% prs.filtwidth : Width of gaussian filter (in samples) to smooth spike train. 
%                 only used for computing % variance explained
% prs.linkfunc  : Name of the link function g ('log','identity' or 'logit').
% prs.lambda    : 1 x N cell array of hyper-parameters for imposing smoothness prior on tuning functions. 
%                 use 0 to impose no prior
% prs.alpha     : Significance level for comparing likelihood values. 
%                 used for model selection
%
% OUTPUT:
% models is a structure containing the results of fitting different model
% variants. It contains the following fields.
% models.class              : M x 1 cell array containing the class of each of the M model variants.
%                             Each cell is a 1 x N logical array that indicates which of the N input
%                             variables were included in the that model variant
% models.testFit            : M x 1 cell array containing results of testing each model variant on the test sets.
%                             Each cell is a k x 6 matrix where the k-th row contains the results of the k-th test set.
%                             The 6 columns contain fraction of variance explained, pearson correlation between actual 
%                             and predicted response, log likelihood, mean-squared prediction error, total number of
%                             spikes, and total number of observations in that order
% models.trainFit           : M x 1 cell array containing results of testing each model variant on the training sets.
% models.wts                : M x 1 cell array of best-fit weights f_i obtained for each of the M model variants.
%                             f_i is empty if x_i was not included in the model variant
% models.x                  : 1 x N cell array of bin centres at which tunings are computed for the input variables.
% models.bestmodel          : The model variant with the highest validated log likelihood.
% models.marginaltunings    : M x 1 cell array of tuning functions for each of the M model variants.
%                             tuning to x_i is empty if x_i was not included in the model variant

%% number of input variables
nvars = length(xt);

%% load analysis parameters
prs = struct2cell(prs);
[xname, xtype, nbins, binrange, lambda, nfolds, dt, filtwidth, linkfunc, alpha] = deal(prs{:});

%% define undefined analysis parameters
if isempty(alpha),      alpha     = 0.05;   end
if isempty(linkfunc),   linkfunc  = 'log';  end
if isempty(filtwidth),  filtwidth = 3;      end
if isempty(nfolds),     nfolds    = 10;     end
if isempty(nbins), nbins = cell(1,nvars); nbins(:) = {10}; end
if isempty(binrange),   binrange  = [];     end
if isempty(dt),         dt        = 0.02;   end

%% compute inverse-link function
if strcmp(linkfunc,'log');          invlinkfunc = @(x) exp(x);
elseif strcmp(linkfunc,'identity'); invlinkfunc = @(x) x;
elseif strcmp(linkfunc,'logit');    invlinkfunc = @(x) exp(x)./(1 + exp(x));
end

%% define bin range
if isempty(binrange), binrange = mat2cell([min(cell2mat(xt));max(cell2mat(xt))],2,strcmp(xtype,'2D')+1); end
% express binrange in units of dt for event type
indx = find(strcmp(xtype,'event'));
for i = indx; binrange{i} = round(binrange{i}/dt); end

%% encode variables in 1-hot format
x    = cell(1,nvars);  % 1-hot representation of xt
xc   = cell(1,nvars);  % bin centres
nprs = cell(1,nvars);  % number of parameters (weights) / bin combinations
% Addition by Jean: observed distributions
observedFR  = cell(1,nvars);
nDataPoints = cell(1,nvars);
for i = 1:nvars
    if size(xt{i},2) > 3
        error('xt wrong dimension...');
    end
    [x{i}, xc{i}, nprs{i}] = Encode1hot(xt{i}, xtype{i}, binrange{i}, nbins{i},dt);
    % special grids for the variables below
    if strcmp(xname{i},'FacingLocation') || strcmp(xname{i},'SpatialView')
        [x{i}, xc{i}, nprs{i}]  = SpatialViewTo1hotConn(xt{i});
%     elseif strcmp(xname{i},'OriSpeed')
%         [x{i}, xc{i}, nprs{i}]  = OriSpeedTo1hot(xt{i});
%     elseif strcmp(xname{i},'EgoBoundary')
%         [x{i}, xc{i}, nprs{i}]  = EBTo1hot(xt{i});
    end
    nDataPoints{i} = sum(x{i});
    observedFR{i}  = sum(x{i}.*repmat(yt,1,size(x{i},2)))./sum(x{i})/dt;
end

models.observedFR  = observedFR;
models.nDataPoints = nDataPoints;
Px = cellfun(@(x) sum(x)/sum(x(:)),x,'UniformOutput',false); % probability of being each state (used for marginalization in the end)

%% define model combinations to fit
nModels    = sum(arrayfun(@(k) nchoosek(nvars,k), 1:nvars));
Xname      = cell(nModels,1);
X          = cell(nModels,1);
Xtype      = cell(nModels,1);
Nprs       = cell(nModels,1);
Lambda     = cell(nModels,1);
ModelCombo = arrayfun(@(k) mat2cell(nchoosek(1:nvars,k),ones(nchoosek(nvars,k),1)),1:nvars,'UniformOutput',false);
ModelCombo = vertcat(ModelCombo{:});
Model      = cell(nModels,1);
for i = 1:nModels; Model{i} = false(1,nvars); end
for i = 1:nModels
    Model{i}(ModelCombo{i}) = true;
    Xname{i}  = xname(Model{i});
%     X{i}      = cell2mat(x(Model{i}));  % X{i} stores inputs for the i^th model, takes lots of memory, no need to preload
    Xtype{i}  = xtype(Model{i});
    Nprs{i}   = cell2mat(nprs(Model{i}));
    Lambda{i} = lambda(Model{i});
end

%% define filter to smooth the firing rate
t = linspace(-2*filtwidth,2*filtwidth,4*filtwidth + 1);
h = exp(-t.^2/(2*filtwidth^2));
h = h/sum(h);

%% fit all models
fprintf(['...... Fitting ' linkfunc '-link model\n']);
models.class    = Model;
models.testFit  = cell(nModels,1);
models.trainFit = cell(nModels,1);
models.wts      = cell(nModels,1);

% Economic version: Instead of fitting model, only initialize NaNs
% fitting in the next section
for i = 1:nModels    
    models.testFit{i}  = nan(nfolds,6);
    models.trainFit{i} = nan(nfolds,6);
    models.wts{i}      = nan(1,sum(Nprs{i})) ;
end

% for n = 1:nModels
%     fprintf('\t- Fitting model %d of %d\n', n, nModels);
%     [models.testFit{n},models.trainFit{n},models.wts{n}] = FitModel(X{n},Xtype{n},Nprs{n},yt,dt,h,nfolds,Lambda{n},linkfunc,invlinkfunc);
% end

models.x = xc;

%% select best model
fprintf('...... Performing forward model selection\n');
testFit = cell2mat(models.testFit);
nrows   = size(testFit,1);
LLvals  = reshape(testFit(:,3), nfolds, nrows/nfolds);  % 3rd column contains likelihood values
% models.bestmodel = ForwardSelect(Model,LLvals,alpha);

ModelMatrix = cell2mat(Model)';
for i = 1:nvars
    if i == 1  % select the best 1-variable model
        modelinx = find(sum(ModelMatrix)==i);
        for n = modelinx
            fprintf('\t- Fitting model %d of %d (%s)\n', n, nModels, char('0'+ModelMatrix(:,n)));
            x_this = cell2mat(x(Model{n}));
            [models.testFit{n},models.trainFit{n},models.wts{n}] = FitModel(x_this,Xtype{n},Nprs{n},yt,dt,h,nfolds,Lambda{n},linkfunc,invlinkfunc,Xname{n}); % added Xname{n} 20200405
        end
        testFit = cell2mat(models.testFit);
        nrows = size(testFit,1);
        LLvals = reshape(testFit(:,3),nfolds,nrows/nfolds);

        bestmodel = find(nanmean(LLvals) == max(nanmean(LLvals(:,modelinx))));
        if isempty(bestmodel)
            break; 
        end
    else  % select the best i-variables model from among those containing the best i-1 variables
        indx1 = (sum(ModelMatrix)==i); % all models containing i variables
        indx2 = (sum(ModelMatrix(ModelMatrix(:,bestmodel)>0,:),1)==i-1);  % all models containing best i-1 variables
        modelinx = find(indx1&indx2);
        for n = modelinx
            fprintf('\t- Fitting model %d of %d (%s)\n', n, nModels,char('0'+ModelMatrix(:,n)));
            x_this = cell2mat(x(Model{n}));
            [models.testFit{n},models.trainFit{n},models.wts{n}] = FitModel(x_this,Xtype{n},Nprs{n},yt,dt,h,nfolds,Lambda{n},linkfunc,invlinkfunc,Xname{n}); % added Xname{n} 20200405
        end
        testFit = cell2mat(models.testFit);
        nrows = size(testFit,1);
        LLvals = reshape(testFit(:,3),nfolds,nrows/nfolds);
                
        bestcandidate = (nanmean(LLvals) == max(nanmean(LLvals(:,modelinx))));
        % significance test :: best new candidate vs current best model
        [pval,~]  = signrank(LLvals(:,bestcandidate),LLvals(:,bestmodel),'tail','right');   % n = 5 folders is at the low end
        [~,pval2] = ttest(LLvals(:,bestcandidate),LLvals(:,bestmodel),'tail','right');
        if pval > alpha && pval2 > alpha, break;
        else, bestmodel = find(bestcandidate); end
    end     
end

if isempty(bestmodel)
    models.bestmodel = NaN;
else
    models.bestmodel = bestmodel;
end

% if isnan(LLvals(1,end))
%     n = size(LLvals,2) ;
%     fprintf('\t- Fitting model %d of %d (%s)\n', n, nModels,char('0'+ModelMatrix(:,n)));
%     x_this = cell2mat(x(Model{n}));
%     [models.testFit{n},models.trainFit{n},models.wts{n}] = FitModel(x_this,Xtype{n},Nprs{n},yt,dt,h,nfolds,Lambda{n},linkfunc,invlinkfunc,Xname{n});
% end

%% match weights 'wts' to corresponding inputs 'x'
models.wts = cellfun(@(x,y) mat2cell(x,1,cell2mat(nprs).*y),models.wts,models.class,'UniformOutput',false);

%% convert weights to response rate (tuning curves) & wrap 2D tunings if any
for i=1:nModels
    for j=1:nvars
        other_vars = (j~=1:nvars & models.class{i});
        if models.class{i}(j)
            if isempty(models.wts{i}(other_vars)) || (strcmp(xtype{j}, 'event')&&all(strcmp(xtype(other_vars),'event'))) % events don't overlap->need not marginalise
                other_factors = 0;
            else
                other_factors = sum(cellfun(@(x,y) sum(x.*y), models.wts{i}(other_vars), Px(other_vars)));
            end
            models.marginaltunings{i}{j} = invlinkfunc(models.wts{i}{j} + other_factors)/dt;
            if strcmp(xtype{j},'2D') && ~strcmp(xname{j},'FacingLocation') && ~strcmp(xname{j},'SpatialView')
                models.marginaltunings{i}{j} = reshape(models.marginaltunings{i}{j},nbins{j}(2),nbins{j}(1));
            end
        else
            models.marginaltunings{i}{j} = [];
        end
    end
end

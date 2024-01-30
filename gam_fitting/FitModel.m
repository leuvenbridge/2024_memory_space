function [testFit,trainFit,param_mean] = FitModel(X,Xtype,Nprs,y,dt,h,nfolds,Lambda,linkfunc,invlinkfunc,Xname)

%% Description
% This function will section the data into nfolds different portions. Each 
% portion is drawn from chunks across the entire recording session. The 
% chunking is done to ensure that each portion includes data from different 
% time points in the experiment (rather than restricted to beginning or the
% end). It will then fit the model to nfolds-1 sections, and test the model 
% performance on the remaining section. This procedure will be repeated 
% nfolds times, with all possible unique testing sections. The fraction of 
% variance explained, the mean-squared error, the log-likelihood increase, 
% and the mean square error will be computed for each test data set. In 
% addition, the learned parameters will be saved for each section of the data.

%% Section the data for k-fold cross-validation and initialise output
nchunks   = 3; % divide data into these many chunks (5 is large enough for experiments that last < 3 hrs)
[~, nprs] = size(X); % number of parameters to fit
nsections = nfolds*nchunks;

% divide the data up into nfolds*nchunks pieces
edges = round(linspace(1,numel(y)+1,nsections+1));

% initialize outputs
testFit  = nan(nfolds,6); % to hold 6 values: var ex, correlation, llh increase, mse, # of spikes, length of test data
trainFit = nan(nfolds,6); % var ex, correlation, llh increase, mse, # of spikes, length of train data
paramMat = nan(nfolds,nprs);

[C,~,ic] = unique(X, 'rows', 'first');
% index of occurrences of each unique state
idx_c = arrayfun(@(i) sort(find(ic==i),'ascend'),1:size(C,1),'UniformOutput',false);

%% perform k-fold cross validation
for k = 1:nfolds
%     fprintf('\t\t- Cross validation fold %d of %d\n', k, nfolds);
    
    %% get test data for the kth fold - comes from chunks across entire session
    test_ind = cell2mat(arrayfun(@(j) edges((j-1)*nfolds + k):edges((j-1)*nfolds + k + 1)-1, 1:nchunks,'UniformOutput',false));
    
    %% get the test_ind from each unique state such that train and test
    % datasets have similar distributions
%     test_ind = [];
%     for istate = 1:length(idx_c)
%         if length(idx_c{istate}) >= nsections
%             edge_this = round(linspace(1, length(idx_c{istate}), nsections+1));
%             idx_ok    = cell2mat(arrayfun(@(j) edge_this((j-1)*nfolds+k):edge_this((j-1)*nfolds+k+1)-1, 1:nchunks, 'UniformOutput', false));
%             idx_ok    = unique(idx_ok);
%             test_ind  = cat(1, test_ind, idx_c{istate}(idx_ok));
%         elseif length(idx_c{istate}) >= k
%             test_ind_this = idx_c{istate}(k:nfolds:end);
%             test_ind  = cat(1, test_ind, test_ind_this);            
%         end
%     end

    test_ind = sort(test_ind, 'ascend');
    
    %%
    test_spikes        = y(test_ind); % test spiking
    smooth_spikes_test = conv(test_spikes,h,'same'); %returns vector same size as original
    smooth_fr_test     = smooth_spikes_test./dt;
    test_X             = X(test_ind,:);
    
    % get training data
    train_ind    = setdiff(1:numel(y),test_ind);
    train_spikes = y(train_ind);
    smooth_spikes_train = conv(train_spikes,h,'same'); %returns vector same size as original
    smooth_fr_train     = smooth_spikes_train./dt;
    train_X = X(train_ind,:);    
    data{1} = train_X;
    data{2} = train_spikes;

    % train the model    
    if strcmp(linkfunc,'log')
        opts = optimset('Algorithm','trust-region','Gradobj','on','Hessian','on','Display','off');
        if k == 1, init_param = 1e-3*randn(nprs, 1); % initialise random parameters for the first training set
        else, init_param = param; end % use final parameters from previous training set
        param = fminunc(@(param) log_link(param,data,Xtype,Nprs,Lambda,Xname),init_param,opts); % fit parameters of log-link model
    elseif strcmp(linkfunc,'identity')
        opts = optimoptions(@fmincon,'Gradobj','on','Hessian','on','Display','off','Algorithm','trust-region-reflective');
        if k == 1, init_param = rand(nprs, 1); % initialise random parameters for the first training set
        else, init_param = param; end % use final parameters from previous training set
        param = fmincon(@(param) identity_link(param,data,Xtype,Nprs,Lambda,Xname),init_param,[],[],[],[],zeros(nprs,1),[],[],opts); % fit parameters of identity-link model
    elseif strcmp(linkfunc,'logit')
        opts = optimset('GradObj','on','Hessian','on','Display','off');
        if k == 1, init_param = 1e-3*randn(nprs, 1); % initialise random parameters for the first training set
        else, init_param = param; end % use final parameters from previous training set
        param = fminunc(@(param) logit_link(param,data,Xtype,Nprs,Lambda,Xname),init_param,opts); % fit parameters of logit-link model
    end
    
    %% %%%%%%%%%%% TEST DATA %%%%%%%%%%%%%
    % compute model predicted firing rate
    fr_hat_test = invlinkfunc(test_X * param)/dt;
    nan_ok      = sum(test_X,2)==0;
    fr_hat_test(nan_ok) = nan;
    smooth_fr_hat_test = conv(fr_hat_test,h,'same'); %returns vector same size as original
    
    % variance explained
    sse = nansum((smooth_fr_hat_test-smooth_fr_test).^2);
    sst = nansum((smooth_fr_test-mean(smooth_fr_test)).^2);
    varExplain_test = 1-(sse/sst);
    
    % linear correlation
    correlation_test = corr(smooth_fr_test,smooth_fr_hat_test,'type','Pearson','rows','complete');
    
    % log-likelihood increase from "mean firing rate model" - NO SMOOTHING
    r_test = invlinkfunc(test_X * param);
    r_test(nan_ok) = [];
    test_spikes(nan_ok) = [];
    n = test_spikes;
    meanFR_test = nanmean(test_spikes);
    log_llh_test_model = nansum(r_test-n.*log(r_test)+log(factorial(n)))/nansum(n); % note: log(gamma(n+1)) will be unstable if n is large (which it isn't here)
    log_llh_test_mean = nansum(meanFR_test-n.*log(meanFR_test)+log(factorial(n)))/nansum(n);
    log_llh_test = (-log_llh_test_model + log_llh_test_mean); % nats/spike
    log_llh_test = log_llh_test/log(2); % convert to bits/spike

    % mean-squared-error
    mse_test = nanmean((smooth_fr_hat_test-smooth_fr_test).^2);
    
    % fill in all the relevant values for the test data from the kth fold
    testFit(k,:) = [varExplain_test correlation_test log_llh_test mse_test nansum(n) numel(test_ind)];
    
    %% %%%%%%%%%%% TRAINING DATA %%%%%%%%%%%
    % compute the smooth firing rate
    fr_hat_train = invlinkfunc(train_X * param)/dt;
    nan_ok       = sum(train_X,2)==0;
    fr_hat_train(nan_ok) = nan;
    smooth_fr_hat_train = conv(fr_hat_train,h,'same'); %returns vector same size as original
    
    % variance explained
    sse = nansum((smooth_fr_hat_train-smooth_fr_train).^2);
    sst = nansum((smooth_fr_train-mean(smooth_fr_train)).^2);
    varExplain_train = 1-(sse/sst);
    
    % correlation
    correlation_train = corr(smooth_fr_train,smooth_fr_hat_train,'type','Pearson','rows','complete');
    
    % log-likelihood
    r_train = invlinkfunc(train_X * param);
    r_train(nan_ok) = [];
    train_spikes(nan_ok) = [];
    n_train = train_spikes;
    meanFR_train = nanmean(train_spikes);   
    log_llh_train_model = nansum(r_train-n_train.*log(r_train)+log(gamma(n_train+1)))/nansum(n_train);
    log_llh_train_mean = nansum(meanFR_train-n_train.*log(meanFR_train)+log(gamma(n_train+1)))/nansum(n_train);
    log_llh_train = (-log_llh_train_model + log_llh_train_mean);
    log_llh_train = log_llh_train/log(2); % convert to bits/spike
    
    % mean-squared-error
    mse_train = nanmean((smooth_fr_hat_train-smooth_fr_train).^2);
    
    % fill in all the relevant values for the training data from the kth fold
    trainFit(k,:) = [varExplain_train correlation_train log_llh_train mse_train nansum(n_train) numel(train_ind)];

    % save the parameters
    paramMat(k,:) = param;
end

param_mean = nanmean(paramMat);

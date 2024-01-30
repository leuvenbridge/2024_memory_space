function [models, bestmodel, VAF, pVAF]  = PlotGAM(models,prs,tit)

if nargin<3, tit='';end
%% Description
% This function will generate three plots:
% 1) log likelihood ratios of each model variant (with standard errors),
% the ratios being taken with respect to a one-parameter null model (constant
% firing rate with no tuning).
% 2) Fraction of variance in neural response explained by each model variant.
% 3) Marginal tuning functions of the best model.

%%
fprintf('...... Plotting results\n');

%% load analysis parameters
prs = struct2cell(prs);
[varname,vartype,~, ~,~, nfolds,~,~,~,~] = deal(prs{:});
nvars = length(varname);

% give each combination of variables a name
nModels = length(models.class);
varlabel = cell(1,nvars); modellabel = cell(1,nModels);
for i=1:nvars
    if strcmp(vartype{i},'2D'), varlabel{i} = varname{i}{1}(1); % use first letter of the variable name to label
    else, varlabel{i} = varname{i}(1); end
end
for i=1:nModels, modellabel{i} = cell2mat(varlabel(models.class{i})); end

%% load model info
testFit = cell2mat(models.testFit);
trainFit = cell2mat(models.trainFit);
nrows = size(testFit,1);
bestmodel = models.bestmodel;
LLvals = reshape(testFit(:,3),nfolds,nrows/nfolds); % 3rd column contains likelihood values
Vexp = reshape(testFit(:,1),nfolds,nrows/nfolds); % 1st column contains variance explained
xvals = models.x;
if ~isnan(bestmodel), fvals = models.marginaltunings{bestmodel}; end

VAF = mean(reshape(trainFit(:,1),nfolds,nrows/nfolds));
pVAF = (VAF(nModels)-VAF((nModels-nvars):(nModels-1)))./(1-VAF((nModels-nvars):(nModels-1))) ;
pVAF=pVAF(end:-1:1) ;

modellabel{1} = 'Position-xy';
modellabel{2} = 'Elevation-z';
modellabel{3} = 'Speed';
modellabel{4} = 'HeadDir';
modellabel{5} = 'AngVel(°/s)';
modellabel{6} = 'PitchAngle(°)';

% plot
figure(1); clf
fs=9;
set(gcf,'Position',[400 500 960 730],'PaperPositionMode','auto','Color','w','InvertHardCopy','off');

axes('Units','Pixels','Position',[0 770 900 20],'FontSize',fs,'FontWeight','bold');hold on
text(0.5,0.5,tit,'FontSize',fs+3,'FontWeight','bold','HorizontalAlignment','Center')
axis([0 1 0 1]);axis off

Nc = 4; % plot N x 4 panels
Nr = 1 + 1 + ceil(nvars/Nc); % plot log-likelihood , var explained , tuning to each variable

% likelihoods
% subplot(Nr,Nc,1:Nc); hold on;
axes('Units','Pixels','Position',[70 580 700 100],'FontSize',fs,'FontWeight','bold');hold on
errorbar(nanmean(LLvals),nanstd(LLvals)/sqrt(nfolds),'ok','linewidth',1);
if not(isnan(bestmodel))
    plot(bestmodel,mean(LLvals(:,bestmodel)),'.r','markersize',20);
end
set(gca,'XLim',[0.5 nModels+0.5]); set(gca,'XTick',1:nModels);
set(gca,'XTickLabel',[]);
% legend('Model performance','Selected model','Null model');
ylabel({'Log likelihood','ratio (bits/spike)'});
h = gca;h.YLim(1)=-0.001;y=h.YLim(2);h.XGrid='on';h.XDir='reverse';


% variance explained
axes('Units','Pixels','Position',[70 440 700 100],'FontSize',fs,'FontWeight','bold');hold on
errorbar(nanmean(Vexp),nanstd(Vexp)/sqrt(nfolds),'ok','linewidth',1);
if not(isnan(bestmodel))
    plot(bestmodel,mean(Vexp(:,bestmodel)),'.r','markersize',20);
end
set(gca,'XLim',[0.5 nModels+.5]); set(gca,'XTick',1:nModels);set(gca,'XTickLabel',[]);

ylabel('Vexp');
h = gca;
if h.YLim(2)>0,
    h.YLim(1)=0;y=h.YLim(1) ;
else
    h.YLim(2)=0;y=h.YLim(1) ;
end
h.XGrid='on';h.XDir='reverse';

for i = 1:nModels
    text(i,y-0.05*abs(h.YLim(2)-h.YLim(1)),modellabel{i},'FontSize',fs,'FontWeight','bold','Rotation',90,'HorizontalAlignment','right','VerticalAlignment','middle')
end

axes('Units','Pixels','Position',[830 440 700*nvars/nModels 240],'FontSize',fs,'FontWeight','bold');hold on
set(bar(pVAF),'FaceColor','k') ;
set(gca,'XLim',[0.5 nvars+0.5]); set(gca,'XTick',1:nvars);
set(gca,'XTickLabel',[]);
ylabel({'Partial correlation'});
set(gca,'YLim',[0 1])
h = gca;h.YLim(1)=0;y=h.YLim(2) ;h.XGrid='on';h.XDir='reverse';
for i = 1:nvars
    text(i,y*-0.05,modellabel{i},'FontSize',fs,'FontWeight','bold','Rotation',90,'HorizontalAlignment','right','VerticalAlignment','middle')
end

% plot tuning functions if the best model is better than the null model
if ~isnan(bestmodel)
    tmp = [] ; 
    for i = 1:nvars, tmp = cat(1,tmp,fvals{i}(:));end
    for i = 1:nvars, tmp = cat(1,tmp,models.observedFR{i}(:));end
else
    tmp = [] ;
    for i = 1:nvars, tmp = cat(1,tmp,models.observedFR{i}(:));end
end
ma = ceil(max(tmp)) ;
mi = floor(min(tmp)) ;
if ~isnan(bestmodel)
    for i=1:nvars
        h = axes('Units','Pixels','Position',[70+(i-1)*180 240 120 120],'FontSize',fs,'FontWeight','bold');hold on
        
        if strcmp(vartype{i},'2D') && ~isempty(fvals{i})
            %             subplot(Nr,Nc,2*Nc+i);
            tmp = fvals{i} ; tmp(isnan(reshape(models.observedFR{i}(:),prs{3}{i})))=NaN;
            imagesc(xvals{i}{1},xvals{i}{2},tmp);
            set(gca,'YDir','reverse');  % Dun: 20180601
            xlabel(varname{i}{1}); ylabel(varname{i}{2});
            axis([xvals{i}{1}([1 end]) xvals{i}{2}([1 end])])
            caxis([mi ma]) ;
        elseif ~isempty(fvals{i})
            %             subplot(Nr,Nc,2*Nc+i);
            plot(xvals{i},fvals{i},'Linewidth',2,'Color','k');
            xlabel(varname{i}); ylabel('Firing rate (spk/s)');
            axis([xvals{i}([1 end]) mi ma])
        else
            axis off
        end
    end
end

% if ~isnan(bestmodel)
for i=1:nvars
    h = axes('Units','Pixels','Position',[70+(i-1)*180 60 120 120],'FontSize',fs,'FontWeight','bold');hold on
    
    if strcmp(vartype{i},'2D') %&& ~isempty(fvals{i})
        %             subplot(Nr,Nc,2*Nc+i);
        b=imagesc(xvals{i}{1},xvals{i}{2},reshape(models.observedFR{i}(:),prs{3}{i}));
        % Dun: 20180601
        set(gca,'YDir','reverse');
        %
        
        xlabel(varname{i}{1}); ylabel(varname{i}{2});
        axis([xvals{i}{1}([1 end]) xvals{i}{2}([1 end])])
        caxis([mi ma]) ;
        
    else %~isempty(fvals{i})
        %             subplot(Nr,Nc,2*Nc+i);
        plot(xvals{i},models.observedFR{i}(:),'Linewidth',2,'Color','k');
        xlabel(varname{i}); ylabel('Firing rate (spk/s)');
        axis([xvals{i}([1 end]) mi ma])
    end
end
% end

% load JL_Colormaps.mat ;
% tmp = interp1(1:64,jet_magenta2,1:0.01:64) ;
% tmp(1,:)=1 ;
% colormap(tmp)

VAF = VAF(end:-1:1) ; pVAF=pVAF(end:-1:1);
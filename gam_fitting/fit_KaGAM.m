function [models, prs] = fit_KaGAM(key,tit)
% key = fetch(mHDC.Result*mHDC.SortingIndex*mHDC.Session*mHDC.Experiment & 'neuron_name="H74M__M"' & 'day="2017-11-28"' & 'protocol="HD_Std"')
% key = key(1) ;

sr=50 ;
[spikes,AAG,XY,LFP,time] = DJ_mHDC_Raw_XY_LFP(key,sr);
Az = AAG(:,1) ;
Omega = JL_gradient_period(Az,1/sr,360) ; Omega(abs(Omega)>300)=NaN ; Omega = JLWF(Omega,10) ;

XY(:,1)=XY(:,1)-580;
XY(:,2)=XY(:,2)-440;
XY=XY/140*25 ;


V = JLWF(gradient(XY')'*30,2);
V = sqrt(sum(V.^2,2)) ;

% clf
% plot(XY(:,1),XY(:,2),'Color',[0 0 0]+0.5) ;

LFP(isnan(LFP))=0 ;
[b,a] = butter(2,12/125,'low') ;[d,c]=butter(2,4/125,'high');
LFP = filter(b,a,filter(d,c,LFP(end:-1:1,:))) ;
LFP = filter(b,a,filter(d,c,LFP(end:-1:1,:))) ;

hilb_eeg = hilbert(LFP); % compute hilbert transform
phase = atan2(imag(hilb_eeg),real(hilb_eeg))*180/pi; %inverse tangent (-pi to pi)
phase = interp1((1:length(phase))/250,phase,time,'nearest') ;
phase = mod(phase,360) ;

prs.varname =  {{'Position-x (cm)' , 'Position-y (cm)'} , 'Linear Speed (cm/s)','Angular Speed (°/s)' , 'Head direction (°)' , 'Theta phase (°)'};
prs.vartype = {'2D'  '1D' '1D'  '1Dcirc'  '1Dcirc'};
prs.nbins = {[20 , 20] , 18 , 18 , 18, 18};
prs.binrange = {[-35 -35;35 35], [0 ; 50],[-300 ; 300], [-180;180],[0 360]}
prs.nfolds = 10;
prs.dt = 0.02;
prs.filtwidth = 3;
prs.linkfunc = 'log';
prs.lambda = {10 , 50 , 50 , 50, 50};
prs.alpha = 0.05;


xt = {XY, V, Omega, Az, phase(:)}; % each of the 3 variables must be T x 1 column vector
yt = spikes(:);

%% fit and plot
tic
models = BuildGAM(xt,yt,prs); % fit
toc
save(['C:\Users\Jean\Dropbox\Houston_Matlab_Analysis\KaGAM\' tit '.mat'],'models','prs','key')
% PlotGAM(models,prs); % plot


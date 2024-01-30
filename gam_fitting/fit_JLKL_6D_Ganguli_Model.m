function [models, prs] = fit_JLKL_6D_Ganguli_Model(spikes,Az,XY,LFP)
% Inputs
% spikes: time stamp of spikes
% Az : head direction (azimuth) in DEGREES, sampled at 50Hz
% XY : head position in cm relative to the center of the arena, sampled at 50Hz
% LFP: LFP signal, sampled at 250 Hz
% Arena_Dimension: adjust so that the arena fits tighlty in a square of +/- this parameter in X and Y dimension

sr=50 ;

Omega = JL_gradient_period(Az,1/sr,360) ; Omega(abs(Omega)>300)=NaN ; Omega = JLWF(Omega,10) ;

V = JLWF(gradient(XY')'*sr,2);
V = V(:,1)*cosd(Az)+V(:,2)*sind(Az);

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


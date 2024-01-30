function [x_1hot,xvals,nprs] = Encode1hot(xt,xtype,binrange,nbins,dt)

% x_1hot stores the state each time point occupies
nt = length(xt);
if strcmp(xtype, 'event')
    basistype = 'raisedcosine';
    binedges = linspace(binrange(1),binrange(2),nbins+1);
    xvals = 0.5*(binedges(1:end-1) + binedges(2:end));  % bin centers
    xt(1:1-binedges(1)) = 0; xt(end-binedges(end):end) = 0;
    x_1hot = zeros(nt, prod(nbins));
    basis  = MakeBasis(basistype, nbins, dt, dt*[binedges(1) binedges(end)]);
    for i = 1:size(basis.y,2)
        xt_conv = conv(xt, basis.y(:,i));
        x_1hot(:,i) = xt_conv(1:end-numel(basis.x)+1);
    end
    x_1hot = circshift(x_1hot, round(basis.x(1)/dt));
elseif strcmp(xtype, '2D')
    % set xt values outside binrange to NaN
    ok = xt(:,1)<binrange(1,1) | xt(:,1)>binrange(2,1); xt(ok,1) = nan;
    ok = xt(:,2)<binrange(1,2) | xt(:,2)>binrange(2,2); xt(ok,2) = nan;
    % initialise with zeros
    binedges1 = linspace(binrange(1,1),binrange(2,1),nbins(1)+1);
    xvals{1}  = 0.5*(binedges1(1:end-1) + binedges1(2:end));  % bin centers
    binedges2 = linspace(binrange(1,2),binrange(2,2),nbins(2)+1);
    xvals{2}  = 0.5*(binedges2(1:end-1) + binedges2(2:end));
    x_1hot    = zeros(nt, prod(nbins));
    % identify index of the ith state and set it to 1
    for i = 1:nt
        if ~isnan(xt(i,1)) && ~isnan(xt(i,2))  % NOTE!! for nans, indx1 and indx2 equals to 1!!
            [~, indx1] = min(abs(xt(i,1)-xvals{1}));
            [~, indx2] = min(abs(xt(i,2)-xvals{2}));
            indx = sub2ind([nbins(2) nbins(1)], nbins(2)-indx2+1, indx1);
            x_1hot(i,indx) = 1;
        end
    end
else
    % set xt values outside binrange to NaN
    ok = xt<binrange(1) | xt>binrange(2); xt(ok) = nan;
    % initialise with zeros
    binedges = linspace(binrange(1),binrange(2),nbins+1);
    xvals = 0.5*(binedges(1:end-1) + binedges(2:end));
    x_1hot = zeros(nt,prod(nbins));
    % identify index of the ith state and set it to 1
    for i = 1:nt
        if ~isnan(xt(i))
            [~, indx] = min(abs(xt(i)-xvals));
            x_1hot(i,indx) = 1;
        end
    end
end
nprs = size(x_1hot,2);
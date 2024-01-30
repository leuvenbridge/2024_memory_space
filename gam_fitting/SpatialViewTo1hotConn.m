function [x_1hot, xval, nprs, M] = SpatialViewTo1hotConn(xt)

%% input xt is the nx3 matrix containing the 3D coordinate at each time point

h = 2200;   % height of the arena, in mm
l = 3500;   % side length of the arena, in mm

d = 250;    % spacing
xx  = -l/2:d:l/2;

% Grid for the floor and ceiling
XYZ = [];
xval = cell(1,2);
for i = 1:length(xx)
    for j = 1:length(xx)
        XYZ(end+1,1:3) = [xx(i) xx(j) h];   
        xval{1}(end+1) = xx(i);
        xval{2}(end+1) = xx(i);
    end
end

for i = 1:length(xx)
    for j = 1:length(xx)
        XYZ(end+1,1:3) = [xx(i) xx(j) 0];
        xval{1}(end+1) = xx(i);
        xval{2}(end+1) = xx(i);
    end
end

% Grid for the 4 walls
for i = 1:length(xx)
    for hh = 0:d:h
        XYZ(end+1,1:3) = [xx(i) l/2 hh];
        xval{1}(end+1) = xx(i);
        xval{2}(end+1) = l/2;
    end
end

for i = 1:length(xx)
    for hh = 0:d:h
        XYZ(end+1,1:3) = [xx(i) -l/2 hh];
        xval{1}(end+1) = xx(i);
        xval{2}(end+1) = -l/2;
    end
end

for i = 1:length(xx)
    for hh = 0:d:h
        XYZ(end+1,1:3) = [l/2 xx(i) hh];
        xval{1}(end+1) = l/2;
        xval{2}(end+1) = xx(i);
    end
end

for i = 1:length(xx)
    for hh = 0:d:h
        XYZ(end+1,1:3) = [-l/2 xx(i) hh];
        xval{1}(end+1) = -l/2;
        xval{2}(end+1) = xx(i);
    end
end


%% 1 hot representation of xt in each state
nt = size(xt,1);
x_1hot = zeros(nt,size(XYZ,1));
for i = 1:nt
    if sum(isnan(xt(i,:)))==0
        [~, idx] = sort(pdist2(XYZ,xt(i,:)));   % idx outputs 1 if xt(i,:) is nan!!
        x_1hot(i,idx(1:2)) = 1;
    end
end
nprs = size(x_1hot,2);

%% plot
% subplot(222)
% D = pdist2(XYZ, XYZ) ; % pairwise distance
% histogram(D(:),3000) ;
% xlim([0 800])
% 
% subplot(121)
% cla; hold on
% plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'or') ;
% axis equal

%% visualize connectivity
% d_thres = d+10;
% x = []; y = []; z = [];
% for i = 1:length(XYZ)
%     for j = i+1:length(XYZ)
%         if D(i,j) <= d_thres
%            x(1:2,end+1) = XYZ([i j],1);
%            y(1:2,end+1) = XYZ([i j],2);
%            z(1:2,end+1) = XYZ([i j],3);           
%         end
%     end
% end
% 
% plot3(x,y,z,'-k')

%% Connectivity Matrix for 3D arena space
% D = pdist2(XYZ, XYZ);
% d_thres = d+10;
% M = D*0 ;
% M(D<d_thres) = -1 ; % set -1 between adjacent bins
% % the diagonal elements of M need to be set so that the sum of each line is 0
% n = size(M,1) ;
% M(logical(eye(n))) = 0 ;
% M(logical(eye(n))) = -sum(M) ;
% 
% % Convert it to a sparse matrix
% M = sparse(M) ;

%% Connectivity Matrix for 2D tilt var. which has nbins = [18 9]
% nbins = [18 9];
% XY = zeros(nbins(2),nbins(1));
% M = zeros(nbins(1)*nbins(2),nbins(1)*nbins(2));
% X = 1:nbins(1); Y = 1:nbins(2);
% xx = repmat(Y',1, nbins(1));
% yy = repmat(X,nbins(2),1);
% D = pdist2([xx(:),yy(:)],[xx(:),yy(:)]);
% M = D.*0;
% M(D==1) = -1;
% n = size(M,1);
% M(logical(eye(n))) = 0;
% M(logical(eye(n))) = -sum(M) ;
% 
% % Convert it to a sparse matrix
% M = sparse(M);

end


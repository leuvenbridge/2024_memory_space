function [x_1hot, xval, nprs, M] = SpatialViewTo1hotConn_JL_DM(xt)

% input xt is the nx3 matrix containing the 3D coordinate at each time point

h = 2120;   % height of the arena, in mm
r = 1650;   % radius of the arena, in mm

d = 200;    % spacing
x = -r:d:r;

% Grid for the floor and ceiling
XYZ = [];
xval = cell(1,2);
for i = 1:length(x)
    for j = 1:length(x)
        XYZ(end+1,1:3) = [x(i) x(j) h];   
        xval{1}(end+1) = x(i);
        xval{2}(end+1) = x(i);
    end
end

for i = 1:length(x)
    for j = 1:length(x)
        XYZ(end+1,1:3) = [x(i) x(j) 0];
        xval{1}(end+1) = x(i);
        xval{2}(end+1) = x(i);
    end
end

I = sum(XYZ(:,1:2).^2,2) > (r-d/4)^2 ;
XYZ(I,:)=[];
xval{1}(I)=[];
xval{2}(I)=[];

% Perimeter of the arena: 2*pi*r
% Number of bins along the perimeter of the arena: round(2*pi*r/d)
% Corresponding angle: 
da = 360/round(2*pi*r/d);
for a = 0:da:360
    for hh = rem(h,d)/2:d:h
        XYZ(end+1,1:3) = [r*cosd(a) r*sind(a) hh];
        xval{1}(end+1) = r*cosd(a);
        xval{2}(end+1) = r*sind(a);
    end
end

%% 1 hot representation of xt in each state
nt = size(xt,1);
x_1hot = zeros(nt,size(XYZ,1));
for i = 1:nt
    if sum(isnan(xt(i,:)))==0
        [~, idx] = min(pdist2(XYZ,xt(i,:)));   % idx outputs 1 if xt(i,:) is nan!!
        x_1hot(i,idx) = 1;
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

%% Connectivity Matrix for 2D tile var. which has nbins = [18 6]
% nbins = [18 6];
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


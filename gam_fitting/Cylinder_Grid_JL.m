h = 2120;   % height of the arena, in mm
r = 1650;   % radius of the arena, in mm

d = 300;
x = -r:d:r;
size(x)

% Grid for the floor and ceiling
XYZ = [];
for i = 1:length(x)
    for j = 1:length(x)
        XYZ(end+1,1:3) = [x(i) x(j) h];        
    end
end

for i = 1:length(x)
    for j = 1:length(x)
        XYZ(end+1,1:3)=[x(i) x(j) 0];        
    end
end

I = sum(XYZ(:,1:2).^2,2) > (r-d/4)^2 ;
XYZ(I,:)=[];

% Perimeter of the arena: 2*pi*r
% Number of bins along the perimeter of the arena: round(2*pi*r/d)
% Corresponding angle: 
da = 360/round(2*pi*r/d);
for a = 0:da:360
    for hh = rem(h,d)/2:d:h
        XYZ(end+1,1:3) = [r*cosd(a) r*sind(a) hh] ;
    end
end

subplot(222)
D = pdist2(XYZ,XYZ) ; % pairwise distance
histogram(D(:),3000) ;
xlim([0 800])

subplot(121)
cla; hold on
plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'or') ;
axis equal

x = []; y = []; z = [];
for i = 1:length(XYZ)
    for j = i+1:length(XYZ)
        if D(i,j) <= 320
           x(1:2,end+1) = XYZ([i j],1);
           y(1:2,end+1) = XYZ([i j],2);
           z(1:2,end+1) = XYZ([i j],3);           
        end
    end
end
figure
plot3(x,y,z,'-k')

%% Connectivity Matrix for Ganguli
M = D*0 ;
M(D<320) = -1 ; % set -1 between adjacent bins
% the diagonal elements of M need to be set so that the sum of each line is 0
n = size(M,1) ;
M(logical(eye(n)))=0 ;
M(logical(eye(n)))=-sum(M) ;

% Convert it to a sparse matrix
M = sparse(M) ;


%%


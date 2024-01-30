function Y = JLWF(X,n)

if size(n,2)==2
    Y=JLWF(JLWF(X,n(1))',n(2))' ;
    return
end
if n == 0, Y=X;return;end

if size(X,2)>1
    Y = zeros(size(X)) ;
    for i = 1:size(X,2)
        Y(:,i) = JLWF(X(:,i),n);
    end
    return
end

Y = zeros(size(X)) ;
S = 0 ;
N = 0 ;
imax = size(X,1) ;
% Initialize
for i = 1:n
   if not(isnan(X(i))), S = S+X(i);N = N+1;end 
end

% Start
for i = 1:imax
    % Add next
    if i+n<=imax
        if not(isnan(X(i+n))), S = S+X(i+n);N = N+1;end 
    end   
    % Remove last
    if i-n-1 >0
         if not(isnan(X(i-n-1))), S = S-X(i-n-1);N = N-1;end 
    end
    if N>0, Y(i) = S/N;else Y(i)=NaN;end
end
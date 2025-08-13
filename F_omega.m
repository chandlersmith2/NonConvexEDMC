function X = F_omega(X_in,samples)
% 
    n = size(X_in,1);
    L = n*(n-1)/2;
    I = samples(:,1);
    J = samples(:,2);
    edgeind = I + (J - 1)*n;
    diagind = (1:n + 1:n*n);
    X_diag= diag(X_in);
    Z1 = X_diag(I)+X_diag(J)-2*X_in(edgeind);
    X = zeros(n,n);
    X(edgeind) = -Z1(1:end);
    X(diagind) =  -sum(X,2);
%     X = X./2;

    

%         
%     sz = size(X_in);
%     I = sub2ind(sz,samples(:,1),samples(:,1));
%     J = sub2ind(sz,samples(:,2),samples(:,2));
%     IJ = sub2ind(sz,samples(:,1),samples(:,2));
%     JI = sub2ind(sz,samples(:,2),samples(:,1));
%     ip = X_in(I) + X_in(J) - X_in(IJ)-X_in(JI);
%     X = zeros(sz(1));
%     X(IJ) = -ip;
%     X(JI) = -ip;
%     rowsum = sum(X);
%     for i=1:sz(1)
%         X(i,i) = -rowsum(i);
%     end

return
   
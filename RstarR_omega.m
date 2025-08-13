function RstarR_X = RstarR_omega(X,samples)

n=size(X,1);

I = samples(:,1);
J = samples(:,2);
edgeind = I + (J - 1)*n;
X_diag= diag(X);
Z1 = X_diag(I)+X_diag(J)-2*X(edgeind);
P_omega_D = zeros(n,n);
P_omega_D(edgeind) = Z1(1:end);

JXJ = -1/2*(P_omega_D - 1/n*ones(n,1)*sum(P_omega_D,2)' - 1/n*sum(P_omega_D,1)'*ones(1,n) + 1/n^2*sum(sum(P_omega_D))*ones(n));
P_JXJ = P_omega(JXJ,samples);
RstarR_X = P_JXJ - diag(sum(P_JXJ,1));

return
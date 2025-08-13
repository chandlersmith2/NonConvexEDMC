function X_0 = R_omega(X,samples)


n = size(X,1);
D = ones(n,1)*diag(X)' + diag(X)*ones(1,n) - 2*X;
P_omega_D = P_omega(D,samples);
X_0 = -(1/2)*(P_omega_D - (1/n)*ones(n,1)*sum(P_omega_D,1)-(1/n)*sum(P_omega_D,2)*ones(1,n) + (1/n^2)*sum(sum(P_omega_D))*ones(n));


return      
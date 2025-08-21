function M_omega_X = M_omega(X,samples,p)
% X is an input variable
% samples are the sampled distances
% p is the bernoulli parameter that was used to sample the distances

n = size(X,1);
M_omega_X = RstarR_omega(X,samples)+(p-1)*(1/2)*(1-1/n+2/n^2)*F_omega(X,samples);

return


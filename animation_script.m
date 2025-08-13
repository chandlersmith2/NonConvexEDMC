

data_loc = "swiss";%,"./data/cow.off","swiss","cities"
rate = 0.05;
num_iter = 100;
rel_thresh = 1*10^-7;
if data_loc == "./data/cow.off" %datapoints P and mesh trg
    [P, trg] = ReadOFF(data_loc,'1');
elseif data_loc == "./data/1k.off"
    [P, trg] = ReadOFF(data_loc,'1'); %datapoints P and mesh trg
elseif data_loc == "swiss"
    load('./data/ptswiss.mat');
    P = pt;
elseif data_loc == "cities"
    load('./data/UScities.mat');
    P = spt(:,1:2); %can make the last arg : to get altitude component
else
    disp 'You should write something here to load your dataset!'
end

n = size(P,1);
%dimension of datapoints
d = size(P,2);
%number of squared distances (useful for one of the sampling modes)
L = n*(n-1)/2;
%make sure that the datapoints have zero mean for reconstruction
P = P - sum(P,1)/n; %centering the data

%build the true gram and distance matrices
X = P*P';%gram matrix
D = ones(n,1)*diag(X)'+diag(X)*ones(1,n)-2*X;%distance matrix
for q=1:n
    D(q,q) = 0.0;
end

[samples,Weight] = Construct_Samples(n,n,L,rate,0,1);
X_0 = (1/rate)*R_omega(X,samples);
M_omega_X = M_omega(X,samples,rate);
[X_approx,P_approx,norm_diffs] = distgeo_animator(X_0,samples,d,num_iter,M_omega_X,rel_thresh,X,rate);

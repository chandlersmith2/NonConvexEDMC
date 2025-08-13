
% This determines which algorithm you're using. Set it to
% 1 for M_omega, 2 for F_omega, 3 for R_omega
descent_mode = 3;

%1 for one step hard thresholding, 2 for Riemannian Resampling Initialization
initialization_type = 1;

%set sampling_mode = 0 for Bernoulli sampling of the entries, sampling_mode
%= 1 for a uniform with replacement sampling model
sampling_mode = 0;

%select 1 if you want to print progress at each iteration
disp_progress = 1;

% Select the datasets
data_array = ["./data/1k.off"];%,"./data/cow.off"];%["./data/cow.off","./data/cow.off","swiss"];%,"cities"];

% Choose the sampling rate p
rate_array = [.04,0.02];%[.1,.07,.05,.03,.02,.01];
% rate_array = [0.02,0.015];
% rate_array = [.05];

lsopts1.maxit = 50;
lsopts1.xtol = 1e-8;
lsopts1.gtol = 1e-8; 
lsopts1.ftol = 1e-10; 
lsopts1.alpha  = 1e-3; 
lsopts1.rho  = 1e-4; 
lsopts1.sigma  = 0.1; 
lsopts1.eta  = 0.8; 

%alternating completion parameters
opts1.maxit = 10;
opts1.r = 1.0;
opts1.printenergy = 0;
opts1.printerror = 0;
opts1.tol = rel_thresh;

%if using the Resampling initialization, give an estimate for your
%coherence parameters and set your partition size
coherence = 50;
partition_size = 3;

% Number of times to run each dataset (with a new randomized subset of
% points)
num_trials = 5;

% Max iterations stopping condition for the algorithms
num_iter = 100;

% relative difference threshold 
rel_thresh = 1*10^-10;

% 3d result array
IPM_results = zeros(length(data_array),length(rate_array),num_trials);
RMSE_results = zeros(length(data_array),length(rate_array),num_trials);

ac_ipm_errs = zeros(length(data_array),length(rate_array),num_trials);
% 
% figure
% tiledlayout(length(data_array),length(rate_array),'Padding', 'none', 'TileSpacing', 'compact');

for i=1:length(data_array)
    data_loc = data_array(i);
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
    elseif data_loc == "ellipse"
        n = 1000;
        theta = 2 * pi * rand(n, 1);
        phi = pi * rand(n, 1);
        x = cos(theta) .* sin(phi);
        y = sin(theta) .* sin(phi);
        z = 0.5 * cos(phi); % Highly degenerate along the z-axis
        P_ellipse = [x, y, z];
        % Add two outlier points along the z-axis
        outlier1 = [0, 0, 10];
        outlier2 = [0, 0, -10];
        P = [P_ellipse; outlier1; outlier2];
    else
        disp 'You should write something here to load your dataset!'
        break
    end
    
    %number of datapoints
    n = size(P,1);
    %dimension of datapoints
    d = size(P,2);
    opts1.rank = d;

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

    % figure
    % tiledlayout(length(rate_array),num_trials) 
    % title('t',data_loc)
    %cycle through rates
    % figure
    % sgtitle(['Dataset: ', char(data_loc), ', Lambda: ', num2str(lambda)])
    % tiledlayout(length(rate_array),num_trials) 
    figure
    sgtitle(['Convergence for Figure ', char(data_loc)])
    tiledlayout(length(rate_array), num_trials)

    for j=1:length(rate_array)
        rate = rate_array(j);
        for k=1:num_trials       

            %build the random samples of the matrix
            %samples is the list passed through to the algorithm, Weight is
            %a binary masking matrix that indicates where entries have been
            %sampled. Use spy(Weight) to see where your entries are
            %plotted!
            if initialization_type ~= 2
                [samples,Weight] = Construct_Samples(n,n,rate,sampling_mode,1);
            else
                [samples,~] = Construct_Samples(n,n,L,rate,sampling_mode,0);
            end

            if descent_mode == 1
                if initialization_type == 1
                    %One step hard thresholding
                    sym_samples = samples;
                    % m = length(samples);
                    % R_omega_X = R_omega(X,sym_samples);
                    % X_0 = L/m*R_omega_X;
                    % F_omega_X = F_omega(X,sym_samples);
                    X_0 = (1/rate)*R_omega(X,sym_samples);
                    % X_0 = RstarR_omega(X,sym_samples);
                    M_omega_X = M_omega(X,sym_samples,rate);
            

                elseif initialization_type == 2
                    %Riemannian resampling
                    X_0 = riemannian_resampling(X,partition_size,samples, ...
                        coherence,d);
                    I = samples(:,1);
                    J = samples(:,2);
                    A = sub2ind([n n],I,J);
                    Weight = zeros(n);
                    Weight(A) = 1;
                    for l=1:n
                        %diagonal entries are known i.e. D_ii = 0
                        Weight(l,l)= 1;
                        for p=l+1:n
                          %make the weight matrix symmetric
                          Weight(p,l)= Weight(l,p);
                        end
                    end
                    [I,J] = find(Weight==1);
                    sym_samples = [I,J];

                    F_omega_X = F_omega(X,sym_samples);
                else
                    disp 'Only hard thresholding and Riemannian resampling are currently implemented'
                    break
                end
               
                % [X_approx,P_approx,norm_diffs] = distgeo_rgrad_Fomega(X_0, ...
                %                                 sym_samples,d,num_iter, ...
                %                                 F_omega_X,rel_thresh,X); 
                [~,X_0,output] = alternating_completion(D,Weight,opts1,lsopts1);
                % ac_ipm_errs(i,j,k) = output.ReconError;
                [X_approx,P_approx,norm_diffs] = distgeo_rgrad_Momega(X_0,sym_samples,d,num_iter,M_omega_X,rel_thresh,X,rate);
                % X_trunc = hard_thresh(X_approx,d,1);
                IPM_err = norm(X_approx-X,'fro')/norm(X,'fro');
                IPM_results(i,j,k) = IPM_err; 
                % P_approx = classical_edg(X_approx,d,0);
                
                P_edg = P_approx(:,1:d);
                P_edg = P_edg - mean(P_edg);
                [U,~,V] = svd(P_edg'*P);
                RR = U*V';
                P_edg_rot = P_edg*RR;
                RMSE = sqrt(mean(sum((P_edg_rot - P).^2,2)));
                RMSE_results(i,j,k) = RMSE;

                % nexttile
                % if data_loc == "./data/1k.off" || data_loc == "./data/cow.off"
                %     ViewMesh(P_edg_rot,trg)
                % elseif data_loc == "swiss"
                %     scatter3(P_edg_rot(:,1),P_edg_rot(:,2),P_edg_rot(:,3),'.')
                %     axis equal
                %     axis off
                % end

                if disp_progress == 1
                    disp("M_omega " + data_loc + "rate " + string(rate) + ...
                        " "+string(k)+"th trial IPM relative difference = " ...
                        + string(IPM_err))
                    disp(string(k)+"th trial RMSE = " + string(RMSE))
                end
                
                % norm_diffs = norm_diffs(norm_diffs>0);
                % nexttile
                % semilogy(1:length(norm_diffs),norm_diffs)
                % xlabel("num iter")
                % ylabel("rate = "+string(rate))


                figure(1) % Ensure plotting on the first tiledlayout figure
                nexttile
                semilogy(norm_diffs,'LineWidth',2)
                title(['Rate = ',num2str(rate),', Trial = ',num2str(k)])
                xlabel('Iteration')
                ylabel('Relative Difference')
                grid on

                figure(2) % Ensure plotting on the second tiledlayout figure
                nexttile
                scatter3(P_edg_rot(:,1), P_edg_rot(:,2), P_edg_rot(:,3), 10, 'blue', 'filled');
                hold on;
                scatter3(P(:,1), P(:,2), P(:,3), 10, 'red', 'filled');
                hold off;
                title(['Rate = ', num2str(rate), ', Trial = ', num2str(k)])
                xlabel('X')
                ylabel('Y')
                zlabel('Z')
                legend('P_{approx}', 'P', 'Location', 'best')
                grid on

            elseif descent_mode == 2
                if initialization_type == 1
                    %One step hard thresholding
                    sym_samples = samples;
                    m = length(samples);
                    R_omega_X = R_omega(X,sym_samples);
                    X_0 = L/m*R_omega_X;
                    F_omega_X = F_omega(X,sym_samples);


                elseif initialization_type == 2
                    %Riemannian resampling
                    X_0 = riemannian_resampling(X,partition_size,samples, ...
                        coherence,d);
                    I = samples(:,1);
                    J = samples(:,2);
                    A = sub2ind([n n],I,J);
                    Weight = zeros(n);
                    Weight(A) = 1;
                    for l=1:n
                        %diagonal entries are known i.e. D_ii = 0
                        Weight(l,l)= 1;
                        for p=l+1:n
                          %make the weight matrix symmetric
                          Weight(p,l)= Weight(l,p);
                        end
                    end
                    [I,J] = find(Weight==1);
                    sym_samples = [I,J];

                    F_omega_X = F_omega(X,sym_samples);
                else
                    disp 'Only hard thresholding and Riemannian resampling are currently implemented'
                    break
                end
               
                [X_approx,P_approx,norm_diffs] = distgeo_rgrad_Fomega(X_0, ...
                                                sym_samples,d,num_iter, ...
                                                F_omega_X,rel_thresh,X); 
                IPM_err = norm(X_approx-X,'fro')/norm(X,'fro');
                IPM_results(i,j,k) = IPM_err; 
                
                P_edg = P_approx(:,1:d);
                P_edg = P_edg - mean(P_edg);
                [U,~,V] = svd(P_edg'*P);
                RR = U*V';
                P_edg_rot = P_edg*RR;
                RMSE = sqrt(mean(sum((P_edg_rot - P).^2,2)));
                RMSE_results(i,j,k) = RMSE;

                if disp_progress == 1
                    disp("M_omega " + data_loc + "rate " + string(rate) + ...
                        " "+string(k)+"th trial IPM relative difference = " ...
                        + string(IPM_err))
                    disp(string(k)+"th trial RMSE = " + string(RMSE))
                end
                
                norm_diffs = norm_diffs(norm_diffs>0);
                nexttile
                semilogy(1:length(norm_diffs),norm_diffs)
                xlabel("num iter")
                ylabel("rate = "+string(rate))

            elseif descent_mode == 3
                if initialization_type == 1
                    sym_samples = samples;
                    m = length(samples);
                    R_omega_X = R_omega(X,samples);
                    X_0 = L/m*R_omega_X;
                elseif initialization_type == 2
                    X_0 = riemannian_resampling(X,partition_size,samples,coherence,d);
                    I = samples(:,1);
                    J = samples(:,2);
                    A = sub2ind([n n],I,J);
                    Weight = zeros(n);
                    Weight(A) = 1;
                    for l=1:n
                        %diagonal entries are known i.e. D_ii = 0
                        Weight(l,l)= 1;
                        for p=l+1:n
                          %make the weight matrix symmetric
                          Weight(p,l)= Weight(l,p);
                        end
                    end
                    [I,J] = find(Weight==1);
                    sym_samples = [I,J];
                    R_omega_X = R_omega(X,sym_samples);
                else
                    disp 'Choose either 0 or 1 for initialization'
                    break
                end
                [~,X_0,output] = alternating_completion(D,Weight,opts1,lsopts1);
                [X_approx,P_approx,norm_diffs] = distgeo_conjgrad_Romega(X_0, ...
                                                            sym_samples,d,num_iter, ...
                                                            0.001,1,...
                                                            R_omega_X,rel_thresh,X);    
                % [X_approx,P_approx,norm_diffs] = distgeo_rgrad_Romega(X_0, ...
                %                                             sym_samples,d,num_iter, ...
                %                                             R_omega_X,rel_thresh,X);    

                IPM_err = norm(X_approx-X,'fro')/norm(X,'fro');
                IPM_results(i,j,k) = IPM_err;

                P_edg = P_approx(:,1:d);
                P_edg = P_edg - mean(P_edg);
                [U,~,V] = svd(P_edg'*P);
                RR = U*V';
                P_edg_rot = P_edg*RR;
                RMSE = sqrt(mean(sum((P_edg_rot - P).^2,2)));

                
                if disp_progress == 1
                    disp("R_omega " + data_loc + "rate " + string(rate) + " " ...
                        +string(k)+"th trial IPM relative difference = " ...
                        + string(IPM_err))
                    disp(string(k)+"th trial RMSE = " + string(RMSE))
                end

                norm_diffs = norm_diffs(norm_diffs>0);
                nexttile
                semilogy(1:length(norm_diffs),norm_diffs)
                xlabel("num iter")
                ylabel("rate = "+string(rate))
                title(data_loc+" trial " + string(k))
            else
                disp('Only three algorithms to choose from!')
            end
        end
    end
end

AC_IPMS_20IL_3OL = IPM_results;
ac_alone_ipm_errs_20IL_3OL = ac_ipm_errs;

% 
% k = 300;  % for example
% 
% % Ensure you don't exceed the number of available points
% k = min(k, size(P_head, 1));
% 
% % Get random indices
% idx = randperm(size(P_head, 1), k);
% 
% % Select the subset
% P_subset = P_head(idx, :);
% 
% % Scatter plot the x and y coordinates
% scatter3(P_subset(:,1), P_subset(:,2), P_subset(:,3), 50, 'filled');
% xlabel('x');
% ylabel('y');
% title('Random Subset of Filtered Points');
% axis equal;
% 
% theta = deg2rad(45);
% 
% % Rotation matrix about z-axis
% Rz = [cos(theta), -sin(theta), 0;
%       sin(theta),  cos(theta), 0;
%       0,           0,          1];
% 
% % Apply the rotation
% P_rotated = (Rz * P_subset')'; 
% 
% scatter3(P_rotated(:,1), P_rotated(:,2),P_rotated(:,3), 50, 'filled');

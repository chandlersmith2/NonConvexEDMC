sampling_mode = 0;

%select 1 if you want to print progress at each iteration
disp_progress = 1;

% Select the datasets
data_array = ["./data/1k.off"];%,"./data/cow.off"];%["./data/cow.off","./data/cow.off","swiss"];%,"cities"];

% Choose the sampling rate p
rate_array = [.03,.02,.01];
% rate_array = [0.02,0.015];
% rate_array = [.05];


% Number of times to run each dataset (with a new randomized subset of
% points)
num_trials = 5;

% Max iterations stopping condition for the algorithms
num_iter = 100;

% relative difference threshold 
rel_thresh = 1*10^-10;

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


% 3d result array
IPM_results = zeros(length(data_array),length(rate_array),num_trials);
RMSE_results = zeros(length(data_array),length(rate_array),num_trials);


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

    figure
    tiledlayout(length(rate_array),num_trials) 
    title('t',data_loc)
    %cycle through rates
    figure
    sgtitle(['Dataset: ', char(data_loc)])
    tiledlayout(length(rate_array),num_trials) 
    figure
    sgtitle(['Convergence for Figure ', char(data_loc)])
    tiledlayout(length(rate_array), num_trials)

    for j=1:length(rate_array)
        rate = rate_array(j);
        for k=1:num_trials       
            [samples,Weight] = Construct_Samples(n,n,rate,sampling_mode,1);
            M_omega_X = M_omega(X,samples,rate);
            [~,X_0,output] = alternating_completion(D,Weight,opts1,lsopts1);
            [X_approx,P_approx,norm_diffs] = distgeo_rgrad_Momega(X_0,samples,d,num_iter,M_omega_X,rel_thresh,X,rate);
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

            
        end
    end
end



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

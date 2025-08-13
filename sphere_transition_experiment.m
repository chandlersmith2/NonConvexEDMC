n=100;
d_list = 2:10;
% rate_array = sort(.02:.01:.1,'descend');
oversampling_rate_array = linspace(1,4,15);
num_trials = 5;
success_thresh = 1e-3;

disp_progress = 1;

% relative difference threshold 
rel_thresh = 1*10^-15;
num_iter = 100;

lsopts1.maxit = 20;
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
IPM_results = zeros(length(d_list),length(oversampling_rate_array),num_trials);
RMSE_results = zeros(length(d_list),length(oversampling_rate_array),num_trials);

% Normalize each point to lie on the sphere


for i=1:length(d_list)
    d = d_list(i);
    opts1.rank = d;
    randomPoints = randn(n, d);
    norms = sqrt(sum(randomPoints.^2, 2));
    P = randomPoints ./ norms;
    
    P = P-mean(P,1);
    
    X = P*P';

    D = ones(n,1)*diag(X)'+diag(X)*ones(1,n)-2*X;%distance matrix
    % 
    % figure
    % tiledlayout(length(oversampling_rate_array),num_trials) 
    for j=1:length(oversampling_rate_array)
        oversampling_val = oversampling_rate_array(j);
        rate = (n*d-d*(d-1)/2)/(n*(n-1)/2)*oversampling_val;
        for k=1:num_trials       
            [samples,Weight] = Construct_Samples(n,n,rate,0,1);
            sym_samples = samples;
            M_omega_X = M_omega(X,sym_samples,rate);
            % R_omega_X = R_omega(X,sym_samples);
            
               
                % [X_approx,P_approx,norm_diffs] = distgeo_rgrad_Fomega(X_0, ...
                %                                 sym_samples,d,num_iter, ...
                %                                 F_omega_X,rel_thresh,X); 
            [~,X_0,output] = alternating_completion(D,Weight,opts1,lsopts1);
            [X_approx,P_approx,norm_diffs] = distgeo_rgrad_Momega(X_0,sym_samples,d,num_iter,M_omega_X,rel_thresh,X,rate);
            % [X_approx,P_approx,norm_diffs] = distgeo_rgrad_Romega(X_0,sym_samples,d,num_iter,R_omega_X,rel_thresh,X);
            X_trunc = hard_thresh(X_approx,d,1);
            IPM_err = norm(X_trunc-X,'fro')/norm(X,'fro');
            IPM_results(i,j,k) = IPM_err; 
            
            P_edg = P_approx(:,1:d);
            P_edg = P_edg - mean(P_edg);
            [U,~,V] = svd(P_edg'*P);
            RR = U*V';
            P_edg_rot = P_edg*RR;
            RMSE = sqrt(mean(sum((P_edg_rot - P).^2,2)));
            RMSE_results(i,j,k) = RMSE;

            if disp_progress == 1
                disp("M_omega dim" + string(d) + "os rate " + string(oversampling_val) + ...
                    " "+string(k)+"th trial IPM relative difference = " ...
                    + string(IPM_err))
                % disp(string(k)+"th trial RMSE = " + string(RMSE)) 
            end
            
            % norm_diffs = norm_diffs(norm_diffs>0);
            % nexttile
            % semilogy(1:length(norm_diffs),norm_diffs)
            % xlabel("num iter")
            % ylabel("rate = "+string(rate))
            

        end
    end
end

figure
count_below_threshold = sum(IPM_results < success_thresh, 3); % Count occurrences below 10^-3 across trials
gca = imagesc(oversampling_rate_array(1:end),d_list,count_below_threshold(:,1:end));
    % set(gca, 'LineStyle','none');
grid off
clim([0, num_trials]);
% Set colormap to grayscale
colormap(gray);
% Add a colorbar
colorbar;
xlabel('Oversampling Rate', 'FontWeight', 'bold', 'FontSize', 14);
ylabel('Sphere Dimension', 'FontWeight', 'bold', 'FontSize', 14);

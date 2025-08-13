% [P, trg] = ReadOFF('./data/1k.off','1');
rng(1)
n=100;
d=3;
randomPoints = randn(n, d);
norms = sqrt(sum(randomPoints.^2, 2));
P = randomPoints ./ norms;

P = P-mean(P,1);

X = P*P';

D = ones(n,1)*diag(X)'+diag(X)*ones(1,n)-2*X;%

n = size(P,1);
d = size(P,2);
L = n*(n-1)/2;
num_iter = 5000;
rel_thresh = 1e-12;
success_thresh = 1e-2;
sampling_mode = 0;
num_trials = 900;
gammas = linspace(-2, -0.9, 10); % Create a list of gamma values from 10^-5 to 1
oversampling_rate_array = linspace(1,5,10); % Generate 10 values between 0.1 and 0.01 in reverse order
results = zeros(length(gammas), length(oversampling_rate_array), num_trials);
count_array = zeros(length(gammas), length(oversampling_rate_array));

total_iters = num_trials*size(gammas,2)*size(oversampling_rate_array,2);
iter_num = 0;

for i = 1:length(gammas)
    gamma = 10^gammas(i); % Get the current gamma value
    disp(['Current gamma: ', num2str(gamma)]);
    noise_level = gamma;
    R = -noise_level + 2 * noise_level * rand(size(P)); % Generate random perturbation
    % X = P*P'; % Compute the Gram matrix
    X_perturbed = X + R*R'; % Perturb the Gram matrixnum_trials
    P_perturbed = P + R; % Perturb the original point cloud P
    P_perturbed = P_perturbed - mean(P_perturbed); % Center the perturbed point cloud
    for j=1:length(oversampling_rate_array)
        oversampling_val = oversampling_rate_array(j);
        rate = (n*d-d*(d-1)/2)/(n*(n-1)/2)*oversampling_val;% Get the current rate value
        disp(['Current oversampling rate: ', num2str(oversampling_val)]);
        for trial = 1:num_trials
            % if mod(trial,50)==0
            %     disp(['Trial number: ', num2str(trial)]);
            % end
            [samples,Weight] = Construct_Samples(n,n,rate,sampling_mode,1);
            X_0 = (1/rate)*R_omega(X_perturbed,samples);
            M_omega_X = M_omega(X_perturbed,samples,rate);
            [X_approx,~,~] = distgeo_rgrad_Momega(X_0,samples,d,num_iter,M_omega_X,rel_thresh,X,rate);
            ipm_err = norm(X_approx-X,'fro')/norm(X,'fro');
            results(i, j, trial) = ipm_err;
            iter_num = iter_num+1;
            
        end
        count_array(i,j) = sum(results(i,j,:)<success_thresh);
        disp('Progress: ' +string(iter_num/total_iters*100))
    end
end

figure
count_below_threshold = sum(results < success_thresh, 3); % Count occurrences below 10^-3 across trials
gca = imagesc(gammas,oversampling_rate_array,count_below_threshold');
% gca = imagesc(gammas,oversampling_rate_array,all500trials');

% gca = imagesc(gammas,oversampling_rate_array,count_thresh');

    % set(gca, 'LineStyle','none');
grid off
clim([0, num_trials]);
% clim([0, 500]);
% Set colormap to grayscale
colormap(gray);
% Add a colorbar
colorbar;
ylabel('Oversampling Rate', 'FontWeight', 'bold', 'FontSize', 14);
xlabel('Noise Exponent (log_{10}(Noise))', 'FontWeight', 'bold', 'FontSize', 14);

beep
results_coarse0to100 = results;
% 
% count_below_threshold_51to100 = count_below_threshold;
% results_51to100 = results;
% % save('full_noise_results',"full_results")
% 
% results_0to100 = zeros(40,40,100);
% results_0to100(:,:,1:50) = results_0to50;
% results_0to100(:,:,51:end) = results_51to100;
% 
% figure
% count_below_threshold = sum(results_0to100 < success_thresh, 3); % Count occurrences below 10^-3 across trials
% gca = imagesc(gammas,oversampling_rate_array,count_below_threshold');
% % gca = imagesc(gammas,oversampling_rate_array,all500trials');
% 
% % gca = imagesc(gammas,oversampling_rate_array,count_thresh');
% 
%     % set(gca, 'LineStyle','none');
% grid off
% clim([0, 100]);
% % clim([0, 500]);
% % Set colormap to grayscale
% colormap(gray);
% % Add a colorbar
% colorbar;
% ylabel('Oversampling Rate', 'FontWeight', 'bold', 'FontSize', 14);
% xlabel('Noise Exponent (log_{10}(Noise))', 'FontWeight', 'bold', 'FontSize', 14);

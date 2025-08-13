function [X_l,P_approx,norm_diffs] = distgeo_animator(X_0,samples,r,max_iter,M_omega_X,rel_thresh,X,p)
%INPUTS:
% X_0: Initialization. In the testing script, this is either the one step
% hard thresholding, or the resampled initialization. Note that the hard
% thresholding step is done in this algorithm, and thus does not need to be
% performed in any pre-processing.
%samples: This is an m by 2 array of sampled indices. It needs to be
%symmetric.
%d: Rank of the manifold for the descent procedure
%max_iter: Max number of iterations before the algorithm is terminated
%F_omega_X: This is the visible information from the sampling, and is
%important for the gradient
%rel_thresh: This is another stopping condition, this one is for relative
%Frobenius norm difference between iterates.
%X: This is the ground truth, and is not necessary for the algorithm. If
%one wishes to compute the convergence plots, returned as norm_diffs, one
%can pass the ground truth in as X.

%OUTPUTS:
%X_l: Final gram matrix recovered
%P_approx: MDS is performed on the final X_l to give the final point cloud
%embedding in d dimensions
%norm_diffs: This returns the relative difference in frobenius norm of X_l
%and X at each iterate l

norm_diffs = zeros(max_iter,1);
figure('Position', [100, 100, 1200, 600]);

% Define the number of frames for the animation
numFrames = max_iter;

% Create a subplot for the first plot
subplot(1, 2, 1);
if r == 3
    h1 = plot3(nan, nan, nan);
else
    h1 = plot(nan, nan);
end
P = classical_edg(X, r, 0);
x = P(:, 1);
y = P(:, 2);
if r == 3
    z = P(:, 3);
    axis([min(x) max(x) min(y) max(y) min(z) max(z)]);
else
    axis([min(x) max(x) min(y) max(y)]);
end

% Create a subplot for the second plot
subplot(1, 2, 2);
h2 = semilogy(nan, nan);
axis([1 max_iter rel_thresh 1]);
xlabel('Iteration');
ylabel('Relative Frobenius Norm Difference');
title('Convergence Plot');

% Preallocate an array to store the frames
frames(numFrames) = struct('cdata', [], 'colormap', []);


%initialize the iteration using the hard thresholding operator
[X_l,~,D_l,U_l]=hard_thresh(X_0,r,1);

for l=1:max_iter
    %Construct the Euclidean gradient
    G_l = (M_omega_X - M_omega(X_l,samples,p));
    %Project the Euclidean gradient onto the tangent space at X_l to
    %produce the Riemannian gradient
    PU_Gl = U_l*(U_l'*G_l);
    Pt_Gl = PU_Gl + PU_Gl' - (PU_Gl*U_l)*U_l';

    %Perform a line search in the tangent space to produce the optimal step
    %size
    Momega_Pt_Gl = M_omega(Pt_Gl,samples,p);
    alpha_l=(norm(Pt_Gl,'fro')^2)/(sum(sum(Pt_Gl.*Momega_Pt_Gl)));
    % disp('step size ' + string(alpha_l))
    
    %This step is a way to more efficiently compute the next iteration
    %following (Wei et al., 2020)
    Z = alpha_l* G_l;
    ZU = Z*U_l;
    Y1 = ZU-(U_l*(U_l'*ZU));
    [Q,R] = qr(Y1,'econ');
    M_l = [D_l + U_l'*ZU,               R' ;
          R            ,            zeros(r)];
    % this amounts to hard thresholding back to the manifold
    [U,D_l] = eig(M_l);
    D_l = diag(D_l);
    D_l = real(D_l);
    [D_l,idx] = sort(D_l,'descend');
    D_l = D_l(1:r);
    U = U(:,idx(1:r));
    D_l = diag(D_l);
    %project back to the SPD cone if there are any negative eigenvalues
    for i=1:r
        if D_l(i,i)<0
            D_l(i,i)=0;
        end
    end

    %build new unitary matrix
    U_l = [U_l Q]*U(:,1:r);
    % Construct new iterate
    X_l1 = X_l;
    X_l = U_l*D_l*U_l';
    dif = norm(X_l-X_l1,'fro')/norm(X_l,'fro');
    disp('Iterate ' + string(l) + ' with relative difference ' + string(dif))
    norm_diffs(l) = dif;

  
    P_approx = classical_edg(X_l,r,0);
    [U,~,V] = svd(P_approx'*P);
    RR = U*V';
    P_approx = P_approx*RR;
    x = P_approx(:,1);
    y = P_approx(:,2);
    if r == 3
        z = P_approx(:,3);
        set(h1, 'XData', x, 'YData', y, 'ZData', z, 'LineStyle', 'none', 'Marker', '.', 'Color', 'k');
    else
        set(h1, 'XData', x, 'YData', y, 'LineStyle', 'none', 'Marker', '.', 'Color', 'k');
    end
    refreshdata(h1);
    drawnow;
    

    % Update the semilog plot with the current norm difference
    set(h2, 'XData', 1:l, 'YData', norm_diffs(1:l));
    refreshdata(h2);
    drawnow;
    frames(l) = getframe(gcf);
    if norm(X_l1-X_l,'fro')/norm(X_l,'fro') < rel_thresh
        numFrames = l;
        break
    end
end
filename = 'animation.gif';
for k = 1:numFrames
    % Convert the frame to an image
    [imind, cm] = rgb2ind(frame2im(frames(k)), 256);
    
    % Write to the .gif file
    if k == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.05);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
    end
end

%Perform MDS to return the point cloud representation




return


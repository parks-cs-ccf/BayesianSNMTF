close all; clear; clc;
%% select the experiment
n_sel_exp = 1;
mc_experiments = {'Experiment-1', ... % baseline experiment
    'Experiment-2', ... % added noisey block
    'Experiment-3'}; % added noise to whole observation matrix (note: in this experiment,
                     % the code provided only adds one noise level instead
                     % of varying levels as in the supp material
mc_sel_exp = mc_experiments{n_sel_exp};

% switch between the experiments
switch mc_sel_exp
    case 'Experiment-1'
        %% generate the data
        %%% dimensions
        % number of boxes
        n_v_box = 4;
        n_h_box = 8;
        
        % number of samples in v box and number of genes in h box
        n_v_per_box = 50;
        n_h_per_box = 100;
        
        % dimensions of X
        n_N = n_v_box*n_v_per_box;
        n_D = n_h_box*n_h_per_box;
        
        % number of cancer types and pathways
        n_K = n_v_box;
        n_R = 4;
        
        %%% initializing matrices
        % X
        % mu values per box
        md_mus = [
            1,  0,  0,  0,  1,  0,  0,  0;
            0, -1,  0,  0,  1,  0,  0, -1;
            0,  0,  0,  1,  0, -1,  0,  0;
            -1,  0,  1,  0,  0,  0,  1,  0;
            ];
        
        % sigma values per box
        md_sigs = [
            2, 1, 1, 1, 2, 1, 1, 1;
            1, 2, 1, 1, 2, 1, 1, 2;
            1, 1, 1, 2, 1, 2, 1, 1;
            2, 1, 2, 1, 1, 1, 2, 1;
            ];
        
        % fill X
        md_X = zeros(n_N, n_D);
        for i = 1:n_v_box
            for j = 1:n_h_box
                % indices for each block
                IDX_v = n_v_per_box*(i-1)+1:n_v_per_box*i;
                IDX_h = n_h_per_box*(j-1)+1:n_h_per_box*j;
                
                md_X(IDX_v, IDX_h) = normrnd(md_mus(i,j),md_sigs(i,j),[n_v_per_box, n_h_per_box]);
            end
        end
        
        % additional noise
        md_sig = 0;
        md_X = md_X + normrnd(0,md_sig,[n_N, n_D]);
        
        % figure; imagesc(md_X); title('X: observation (N X D)'); colorbar; caxis([-6,6]);
        
        % U
        md_U = zeros(n_N, n_K);
        for i = 1:n_K
            md_U((i-1)*n_v_per_box + 1:i*n_v_per_box, i) = 1.0;
        end
        
        %figure; imagesc(md_U); colorbar;
        
        % V0
        md_V0 = zeros(n_D, n_R);
        
        for i = 1:n_h_box
            for j = 1:n_v_box
                % indices for each block
                IDX_v = n_h_per_box*(i-1)+1:n_h_per_box*i;
                IDX_h = j;
                
                md_V0(IDX_v, IDX_h) = md_mus(j,i) ~= 0;
                
            end
        end
        
        %figure; imagesc(md_V0'); colorbar;
        
        % gene-gene interaction
        % md_A = zeros(n_D, n_D);
        % for i = 1:n_h_box
        %     IDX = n_h_per_box*(i-1)+1:n_h_per_box*i;
        %     md_A(IDX, IDX) = 1;
        % end
        
        md_ones = ones(n_D, n_D);
        IDX = find(tril(md_ones, -1));
        rand_IDX = randperm(length(IDX));
        d_rate = 0.1;
        sub_IDX = IDX(rand_IDX(1:ceil(d_rate*length(IDX))));
        md_A = zeros(n_D, n_D);
        md_A(sub_IDX) = 1;
        md_A = md_A + md_A';
        
        %figure; imagesc(md_A), colorbar;
        
        %% Bayesian semi-nonnegative matrix tri-factorization
        mc_opt.MaxIters = 30;
        mc_Solution = myfunc_BSNMTF_ssp_final(md_X, md_U, md_V0, md_A, mc_opt);
        
        md_X_hat = mc_Solution.U * mc_Solution.S * (mc_Solution.Z .* mc_Solution.V)';
        
        BSNTMF_plot(md_X, md_U, md_V0, md_A, md_X_hat, mc_Solution.S, mc_Solution.Z, mc_Solution.V, 1);
        
        %% Nonnegative matrix tri-factorization
        md_X_new = [max(md_X, 0), max(-md_X, 0)];
        md_V0_new = [md_V0; md_V0];
        md_A_new = [md_A, md_A; md_A, md_A];
        
        [mr_bestLV, mr_bestLS] = myfuc_NMTF_best_lambda(md_X_new, md_U, md_V0_new, md_A_new);
        [S_Ntri, V_Ntri] = myfunc_NMTF_lambdas(md_X_new, md_U, md_V0_new, md_A_new, mr_bestLV, mr_bestLS);
        
        % [S_Ntri, V_Ntri] = myfunc_NMTF(md_X_new, md_U, md_V0_new, md_A_new);
        
        md_X_tri = mc_Solution.U * (S_Ntri * V_Ntri');
        
        NMTF_plot(md_X_new, md_U, md_V0_new, md_A_new, md_X_tri, S_Ntri, V_Ntri, 1);
        
    case 'Experiment-2'
        %% generate the data
        %%% dimensions
        % number of boxes
        n_v_box = 4;
        n_h_box = 8;
        
        % number of samples in v box and number of genes in h box
        n_v_per_box = 50;
        n_h_per_box = 100;
        
        % dimensions of X
        n_N = n_v_box*n_v_per_box;
        n_D = n_h_box*n_h_per_box;
        
        % number of cancer types and pathways
        n_K = n_v_box;
        n_R = 4;
        
        %%% initializing matrices
        % X
        % mu values per box
        md_mus = [
            1,  0,  0,  0,  1,  0,  0,  0;
            0, -1,  0,  0,  1,  0,  0, -1;
            0,  0,  0,  1,  0, -1,  0,  0;
            -1,  0,  1,  0,  0,  0,  1,  0;
            ];
        
        % sigma values per box
        md_sigs = [
            2, 5, 1, 1, 2, 1, 1, 1;
            1, 2, 1, 1, 2, 1, 1, 2;
            1, 1, 1, 2, 1, 2, 1, 1;
            2, 1, 2, 1, 1, 1, 2, 1;
            ];
        
        % fill X
        md_X = zeros(n_N, n_D);
        for i = 1:n_v_box
            for j = 1:n_h_box
                % indices for each block
                IDX_v = n_v_per_box*(i-1)+1:n_v_per_box*i;
                IDX_h = n_h_per_box*(j-1)+1:n_h_per_box*j;
                
                md_X(IDX_v, IDX_h) = normrnd(md_mus(i,j),md_sigs(i,j),[n_v_per_box, n_h_per_box]);
            end
        end
        
        % additional noise
        md_sig = 0;
        md_X = md_X + normrnd(0,md_sig,[n_N, n_D]);
        
        % figure; imagesc(md_X); title('X: observation (N X D)'); colorbar; caxis([-6,6]);
        
        % U
        md_U = zeros(n_N, n_K);
        for i = 1:n_K
            md_U((i-1)*n_v_per_box + 1:i*n_v_per_box, i) = 1.0;
        end
        
        %figure; imagesc(md_U); colorbar;
        
        % V0
        md_V0 = zeros(n_D, n_R);
        
        for i = 1:n_h_box
            for j = 1:n_v_box
                % indices for each block
                IDX_v = n_h_per_box*(i-1)+1:n_h_per_box*i;
                IDX_h = j;
                
                md_V0(IDX_v, IDX_h) = md_mus(j,i) ~= 0;
                
            end
        end
        
        %figure; imagesc(md_V0'); colorbar;
        
        % gene-gene interaction
        % md_A = zeros(n_D, n_D);
        % for i = 1:n_h_box
        %     IDX = n_h_per_box*(i-1)+1:n_h_per_box*i;
        %     md_A(IDX, IDX) = 1;
        % end
        
        md_ones = ones(n_D, n_D);
        IDX = find(tril(md_ones, -1));
        rand_IDX = randperm(length(IDX));
        d_rate = 0.1;
        sub_IDX = IDX(rand_IDX(1:ceil(d_rate*length(IDX))));
        md_A = zeros(n_D, n_D);
        md_A(sub_IDX) = 1;
        md_A = md_A + md_A';
        
        %figure; imagesc(md_A), colorbar;
        
        %% Bayesian semi-nonnegative matrix tri-factorization
        mc_opt.MaxIters = 30;
        mc_Solution = myfunc_BSNMTF_ssp_final(md_X, md_U, md_V0, md_A, mc_opt);
        
        md_X_hat = mc_Solution.U * mc_Solution.S * (mc_Solution.Z .* mc_Solution.V)';
        
        BSNTMF_plot(md_X, md_U, md_V0, md_A, md_X_hat, mc_Solution.S, mc_Solution.Z, mc_Solution.V, 2);
        
        
        %% Nonnegative matrix tri-factorization
        md_X_new = [max(md_X, 0), max(-md_X, 0)];
        md_V0_new = [md_V0; md_V0];
        md_A_new = [md_A, md_A; md_A, md_A];
        
        [mr_bestLV, mr_bestLS] = myfuc_NMTF_best_lambda(md_X_new, md_U, md_V0_new, md_A_new);
        [S_Ntri, V_Ntri] = myfunc_NMTF_lambdas(md_X_new, md_U, md_V0_new, md_A_new, mr_bestLV, mr_bestLS);
        
        % [S_Ntri, V_Ntri] = myfunc_NMTF(md_X_new, md_U, md_V0_new, md_A_new);
        
        md_X_tri = mc_Solution.U * (S_Ntri * V_Ntri');
        
        NMTF_plot(md_X_new, md_U, md_V0_new, md_A_new, md_X_tri, S_Ntri, V_Ntri, 2);
        
    case 'Experiment-3'
        %% generate the data
        %%% dimensions
        % number of boxes
        n_v_box = 4;
        n_h_box = 8;
        
        % number of samples in v box and number of genes in h box
        n_v_per_box = 50;
        n_h_per_box = 100;
        
        % dimensions of X
        n_N = n_v_box*n_v_per_box;
        n_D = n_h_box*n_h_per_box;
        
        % number of cancer types and pathways
        n_K = n_v_box;
        n_R = 4;
        
        %%% initializing matrices
        % X
        % mu values per box
        md_mus = [
            1,  0,  0,  0,  1,  0,  0,  0;
            0, -1,  0,  0,  1,  0,  0, -1;
            0,  0,  0,  1,  0, -1,  0,  0;
            -1,  0,  1,  0,  0,  0,  1,  0;
            ];
        
        % sigma values per box
        md_sigs = [
            2, 1, 1, 1, 2, 1, 1, 1;
            1, 2, 1, 1, 2, 1, 1, 2;
            1, 1, 1, 2, 1, 2, 1, 1;
            2, 1, 2, 1, 1, 1, 2, 1;
            ];
        
        % fill X
        md_X = zeros(n_N, n_D);
        for i = 1:n_v_box
            for j = 1:n_h_box
                % indices for each block
                IDX_v = n_v_per_box*(i-1)+1:n_v_per_box*i;
                IDX_h = n_h_per_box*(j-1)+1:n_h_per_box*j;
                
                md_X(IDX_v, IDX_h) = normrnd(md_mus(i,j),md_sigs(i,j),[n_v_per_box, n_h_per_box]);
            end
        end
        
        % additional noise
        md_sig = 5;
        md_X = md_X + normrnd(0,md_sig,[n_N, n_D]);
        
        % figure; imagesc(md_X); title('X: observation (N X D)'); colorbar; caxis([-6,6]);
        
        % U
        md_U = zeros(n_N, n_K);
        for i = 1:n_K
            md_U((i-1)*n_v_per_box + 1:i*n_v_per_box, i) = 1.0;
        end
        
        %figure; imagesc(md_U); colorbar;
        
        % V0
        md_V0 = zeros(n_D, n_R);
        
        for i = 1:n_h_box
            for j = 1:n_v_box
                % indices for each block
                IDX_v = n_h_per_box*(i-1)+1:n_h_per_box*i;
                IDX_h = j;
                
                md_V0(IDX_v, IDX_h) = md_mus(j,i) ~= 0;
                
            end
        end
        
        %figure; imagesc(md_V0'); colorbar;
        
        % gene-gene interaction
        % md_A = zeros(n_D, n_D);
        % for i = 1:n_h_box
        %     IDX = n_h_per_box*(i-1)+1:n_h_per_box*i;
        %     md_A(IDX, IDX) = 1;
        % end
        
        md_ones = ones(n_D, n_D);
        IDX = find(tril(md_ones, -1));
        rand_IDX = randperm(length(IDX));
        d_rate = 0.1;
        sub_IDX = IDX(rand_IDX(1:ceil(d_rate*length(IDX))));
        md_A = zeros(n_D, n_D);
        md_A(sub_IDX) = 1;
        md_A = md_A + md_A';
        
        %figure; imagesc(md_A), colorbar;
        
        %% Bayesian semi-nonnegative matrix tri-factorization
        mc_opt.MaxIters = 30;
        mc_Solution = myfunc_BSNMTF_ssp_final(md_X, md_U, md_V0, md_A, mc_opt);
        
        md_X_hat = mc_Solution.U * mc_Solution.S * (mc_Solution.Z .* mc_Solution.V)';
        
        BSNTMF_plot(md_X, md_U, md_V0, md_A, md_X_hat, mc_Solution.S, mc_Solution.Z, mc_Solution.V, 3);
        
        %% Nonnegative matrix tri-factorization
        md_X_new = [max(md_X, 0), max(-md_X, 0)];
        md_V0_new = [md_V0; md_V0];
        md_A_new = [md_A, md_A; md_A, md_A];
        
        [mr_bestLV, mr_bestLS] = myfuc_NMTF_best_lambda(md_X_new, md_U, md_V0_new, md_A_new);
        [S_Ntri, V_Ntri] = myfunc_NMTF_lambdas(md_X_new, md_U, md_V0_new, md_A_new, mr_bestLV, mr_bestLS);
        
        % [S_Ntri, V_Ntri] = myfunc_NMTF(md_X_new, md_U, md_V0_new, md_A_new);
        
        md_X_tri = mc_Solution.U * (S_Ntri * V_Ntri');
        
        NMTF_plot(md_X_new, md_U, md_V0_new, md_A_new, md_X_tri, S_Ntri, V_Ntri, 3);
        
end

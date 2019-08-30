%----------------------------------------------------------------------
% Bayesian semi-nonnegative matrix tri-factorization
%
%
% This code is an implementation of the experiment on 
% the simulation data in the main paper 
% 
% programed by Sunho Park (parks@ccf.org)
%----------------------------------------------------------------------

clear; close all; clc;

mn_K = 3;
mn_NinK = 100;

mn_D = 400;

mn_N = mn_K*mn_NinK;
mn_R = 3;

%- Construct the input matrix X of R^(NxD) 
mm_X = 0.1*randn(mn_N, mn_D);

% settings for the pattern blocks (please see Figure 2 in the main paper)
mr_mu_pos = 2.0;
mr_mu_neg = -mr_mu_pos;

mr_sig = sqrt(2.0);

% group: 1
mm_RanDN =  mr_mu_pos + mr_sig*randn(100, 100);
mm_X(1:100, 1:100) = mm_RanDN;

% group: 2
mm_RanDP = mr_mu_pos + mr_sig*randn(100, 100);
mm_RanDN = mr_mu_neg + mr_sig*randn(100, 100);

mm_X(100+(1:100), 100+(1:100)) = mm_RanDP;
mm_X(100+(1:100), 300+(1:100)) = mm_RanDN;

% group: 3
mm_RanDN = mr_mu_neg + mr_sig*randn(100, 100);
mm_X(201:300, 200+(1:100)) = mm_RanDN;

% plot 
h = figure; 
imagesc(mm_X); title('X: observation (N X D)'); colorbar; caxis([-6, 6]);


%- U of R(NxK): fixed 
mm_U = zeros(mn_N, mn_K);
for mn_i = 1:mn_K
    mv_idx = ((mn_i-1)*mn_NinK + 1):(mn_i*mn_NinK);
    mm_U(mv_idx, mn_i) = 1;
end

%- Initial pathway informatin: Z^0 \in R(DxR)
mm_Z0 = zeros(mn_D, mn_R);

% r = 1
mm_Z0(1:100, 1) = 1;

% r = 2
mm_Z0(200+(1:100), 2) = 1;

% r = 3
mv_Refidx = find(rand(100, 1) > 0.2);

mv_RefidxC = (1:100);
mv_RefidxC = setdiff(mv_RefidxC, mv_Refidx);

mm_Z0([100+mv_RefidxC , 300+(1:100)], 3) = 1;

%- gene-gene interaction network: A \in R^{DxD)
mm_Ones = ones(mn_D, mn_D);
mv_IDX = find(tril(mm_Ones, -1));

mv_RandIDX = randperm(length(mv_IDX));
mr_Rate = 0.1;
mv_SubIDX = mv_IDX(mv_RandIDX(1:round(mr_Rate*length(mv_IDX))));  

mm_GGINet_mat = zeros(mn_D, mn_D);
mm_GGINet_mat(mv_SubIDX) = 1;

mv_chk = find(sum(mm_GGINet_mat, 2) == 0);
for mn_i = 1:length(mv_chk)
    mn_pos = mv_chk(mn_i);
    
    if mn_pos == 1
        continue; 
    end
    
    mv_randidx = randperm(mn_pos-1);    
    mm_GGINet_mat(mn_pos, mv_randidx(1)) = 1;
end

mm_GGINet_mat = mm_GGINet_mat + mm_GGINet_mat';

%--- Bayesian semi-nonnegative matrix tri-factorization
mc_opt.MaxIters = 30;
mc_Solution = bsnmtf(mm_X, mm_U, mm_Z0, mm_GGINet_mat, mc_opt);

%- generate plots

% initial pathway information
figure;
imagesc(mm_Z0); title('Z0: membership');

% reconstructed matrix
mm_X_hat = mc_Solution.U * mc_Solution.S * (mc_Solution.Z .* mc_Solution.V)';

figure; 
imagesc(mm_X_hat); title('X_hat: reconstructed');

% S_hat
disp(mc_Solution.S)

% Z_hat
figure;
imagesc(mc_Solution.Z'); title('V: Pattenrs');

% V_hat
figure;
imagesc(mc_Solution.V'); title('V: Pattenrs');


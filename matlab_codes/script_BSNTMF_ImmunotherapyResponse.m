%--------------------------------------------------------------------------
% Bayesian semi-nonnegative matrix tri-factorization to identify pathways
% associated with cancer phenotypes
%
% Experiment on the metastatic gastric cancer immunotherapy clinical-trial dataset
% please see Sectioni 4.2 in our main paper
%
% Sunho Park (parks@ccf.org)
%--------------------------------------------------------------------------

clc; clear; close all;

disp('Immunotherapy response data')

mstr_pathway = 'PathDB';
disp(['Pathway DB: ', mstr_pathway]);

%--- load the metastatic gastric cancer immunotherapy resonse data
tmp = load('./data/MGC_Imuno_response.mat');

mc_MGA_pID  = tmp.mc_data.pid;
mc_MGA_pInfo  = tmp.mc_data.PInfo;

%- construct the indicator matrix U (N X 2 = {response vs non-response})
mc_ReponseInfo = mc_MGA_pInfo(:, 9);
mc_str_reps = {'CR', 'PR'}; 
mc_str_noreps = {'SD', 'PD'}; 

mm_U = zeros(length(mc_MGA_pID), 2);
mm_U(ismember(mc_ReponseInfo, mc_str_reps), 1) = 1;
mm_U(ismember(mc_ReponseInfo, mc_str_noreps), 2) = 1;

%- load the gene expression data (the observation matrix X: N X D)
mm_XGE_org = tmp.mc_data.X;
mn_N = length(mc_MGA_pID);

mm_XGE_sub = zscore(log(mm_XGE_org + 1));
mc_gene_GE = tmp.mc_data.genes(:, 2);

%- load the BioGrid gene-gene interaction network (A: D X D)
disp('Gene-Gene Interaction Network: bioGrid');
tmp = load('./data/BioGrid_3_4_153.mat');

mc_GGNet.gene = tmp.mc_BioGrid.GeneSymbols;
mc_GGNet.Mat = tmp.mc_BioGrid.Mat;

[mc_ComGenes, ~, mv_Bidx] = intersect(mc_gene_GE, mc_GGNet.gene);
mm_GGNetMat = mc_GGNet.Mat(mv_Bidx, mv_Bidx);

mv_idx = sum(mm_GGNetMat,1)' > 0 & sum(mm_GGNetMat,2) > 0;
mc_Geneymbols = mc_ComGenes(mv_idx);

mn_D = length(mc_Geneymbols);

%- match the genes in both sets 
mm_XGE = NaN*ones(mn_N, mn_D);
[mc_dummy4, mv_Aidx, mv_Bidx] = intersect(mc_Geneymbols, mc_gene_GE); 
mm_XGE(:,mv_Aidx) = mm_XGE_sub(:,mv_Bidx);

mm_GGINet = zeros(mn_D, mn_D);
[mc_dummy5, mv_Aidx, mv_Bidx] = intersect(mc_Geneymbols, mc_GGNet.gene); 
mm_GGINet(mv_Aidx,mv_Aidx) = mc_GGNet.Mat(mv_Bidx,mv_Bidx);

mm_GGINet = sparse(mm_GGINet);
 
%- Construct the pathway membership matrix (Z^0: D X R)
tmp = load('./data/gene_name_info_hg18');
hg18_gene_name = tmp.hg18_gene_name(:, 1);

tmp = load('./data/bipartite_PPI_module');
module_idx_org = tmp.module_idx;

module_idx = module_idx_org;

mc_PathwayInfo = cell(size(module_idx,1),1);
for mn_i = 1:size(module_idx,1)
    mc_PathwayInfo{mn_i} = sprintf('%4dth pathway', mn_i);
end

mm_V0 = sparse([], [], [], mn_D, size(mc_PathwayInfo,1));
[mc_dummy6, mv_Aidx, mv_Bidx] = intersect(mc_Geneymbols, hg18_gene_name);
mm_V0(mv_Aidx,:) = module_idx(:,mv_Bidx)';

mn_R = size(mm_V0,2);

%--- call BSNMTF
mc_opt.MaxIters = 30;
mc_Ctypes = {'reponse', 'non-response'};

mc_Solution = bsnmtf(mm_XGE, mm_U, mm_V0, mm_GGINet, mc_opt);

save(['MGA_', mstr_pathway, '.mat'], 'mc_Solution', 'mc_PathwayInfo', 'mm_V0', 'mc_Geneymbols', '-v7.3');              


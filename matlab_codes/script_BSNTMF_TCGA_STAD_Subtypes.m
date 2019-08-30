%--------------------------------------------------------------------------
% Bayesian semi-nonnegative matrix tri-factorization to identify pathways
% associated with cancer phenotypes
%
% Experiment on TCGA gastric cancer dataset
% please see Sectioni 4.2 in our main paper
%
% Sunho Park (parks@ccf.org)
%--------------------------------------------------------------------------

clc; clear; close all;

%--- TCGA data loading
disp('TCGA STAD subtypes')

mstr_pathway = 'PathDB';
disp(['Pathway DB: ', mstr_pathway]);

tmp = load('./data/TCGA_Data_Compact_10032017.mat');

% select only the gatric cancer samples
mn_Dtmp = length(tmp.mc_TCGA_Data.RNA_genes);
mv_ChkIDX = strcmpi(tmp.mc_TCGA_Data.CTypes_Vec, 'stad') ...
          & (sum(~isnan(tmp.mc_TCGA_Data.RNA_X),2) > 0.8*mn_Dtmp);

mc_PIDs_org = tmp.mc_TCGA_Data.PIDs(mv_ChkIDX);
%- load the gene expression data (the observation matrix X: N X D)
mm_XGE_org = tmp.mc_TCGA_Data.RNA_X(mv_ChkIDX, :);
mc_gene_GE_org = tmp.mc_TCGA_Data.RNA_genes;

%- construct the indicator matrix U (N X 4 = {'CIN' vs 'EBV' vs 'GS' vs 'MSI'})
%- load the gastric molecular subtype information 
fid = fopen(['./data/data_clinical.txt'], 'r');

tline = fgets(fid);
m_cheaders = strsplit(tline,'\t');
m_nLen = length(m_cheaders);

m_strpattern = [repmat('%s',[1,m_nLen]), '%[^\n\r]'];
m_mData = textscan(fid,m_strpattern,'Delimiter','\t');

fclose(fid);

mc_SubType_org = cell(length(mc_PIDs_org), 1);
[mc_dummy1, mv_Aidx, mv_Bidx] = intersect(mc_PIDs_org, m_mData{1,2});

mc_SubType_org(mv_Aidx) = m_mData{1,3}(mv_Bidx);
mc_SubType = cellfun(@(x) sprintf('%s',x), mc_SubType_org,'UniformOutput',false);

mc_Ctypes = {'CIN', 'EBV', 'GS', 'MSI'};

mm_Uorg = zeros(length(mc_PIDs_org), length(mc_Ctypes));
for mn_i = 1:length(mc_Ctypes)
    mm_Uorg(strcmpi(mc_SubType, mc_Ctypes{mn_i}), mn_i) = 1;
end

mv_selidx = sum(mm_Uorg, 2) > 0;

mc_PIDs = mc_PIDs_org(mv_selidx);
mn_N = length(mc_PIDs);

mm_U = mm_Uorg(mv_selidx, :);
mm_XGE_sub_1 = mm_XGE_org(mv_selidx, :);

% To make sure that there is no missing value in the matrix.
mv_ChkIDX02 = sum(isnan(mm_XGE_sub_1),1)==0;

mm_XGE_sub = mm_XGE_sub_1(:,mv_ChkIDX02);
mc_gene_GE = mc_gene_GE_org(mv_ChkIDX02);

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
hg18_gene_name = tmp.hg18_gene_name(:,1);

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
mc_Solution = bsnmtf(mm_XGE, mm_U, mm_V0, mm_GGINet, mc_opt);

save(['TCGA_STAD_', mstr_pathway, '.mat'], 'mc_Solution', 'mc_PathwayInfo', 'mm_V0', 'mc_Geneymbols', '-v7.3');    

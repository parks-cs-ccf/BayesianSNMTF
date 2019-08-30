clc; clear; close all;

disp('Immunotherapy response data')

mn_PathNum = 4;
mstr_fullFiles = {'kegg', 'reactome', 'biocarta', 'PathDB', 'hallmark'};
mstr_pathway = mstr_fullFiles{mn_PathNum};

disp(['Pathway DB: ', mstr_pathway]);

%--- MGA data loading
mstr_path = 'E:\Research_Mar092017\Code\mRNA_Yonsei\ConsensusClustering\paper_results\Manscript_Tables_Figures\Ref\';

tmp = load([mstr_path, 'MGC_Imuno_response.mat']);

mc_MGA_pID  = tmp.mc_data.pid;
mc_MGA_pInfo  = tmp.mc_data.PInfo;

%- construct a cancer type matrix (samples X cancer types)
mc_ReponseInfo = mc_MGA_pInfo(:, 9);
mc_str_reps = {'CR', 'PR'}; 
mc_str_noreps = {'SD', 'PD'}; 

mm_U = zeros(length(mc_MGA_pID), 2);
mm_U(ismember(mc_ReponseInfo, mc_str_reps), 1) = 1;
mm_U(ismember(mc_ReponseInfo, mc_str_noreps), 2) = 1;

%- Read mRNA (gene expression) data
mm_XGE_org = tmp.mc_data.X;

mn_N = length(mc_MGA_pID);

% To make sure that no missing value is in the matrix.'
mm_XGE_sub = zscore(log(mm_XGE_org + 1));
mc_gene_GE = tmp.mc_data.genes(:, 2);

%- Read BioGrid PPI network
mstr_fullPPINet = {'ppiMat', 'bioGrid'};
mstr_PPINetSel = mstr_fullPPINet{2};

disp(['Gene-Gene Interaction Network: ', mstr_PPINetSel]);
if strcmpi(mstr_PPINetSel, 'ppiMat')    
    tmp = load('E:\BackUp\Code\NonnegativeMatrix\Mutation_05282013\data\ppiMatrixTF.mat');
    
    mc_GGNet.gene = tmp.gene_name_ge;
    mc_GGNet.Mat = tmp.ppiMatrixTF;
else
    tmp = load('./data/BioGrid_3_4_153.mat');
            
    mc_GGNet.gene = tmp.mc_BioGrid.GeneSymbols;
    mc_GGNet.Mat = tmp.mc_BioGrid.Mat;
end

[mc_ComGenes, ~, mv_Bidx] = intersect(mc_gene_GE, mc_GGNet.gene);
mm_GGNetMat = mc_GGNet.Mat(mv_Bidx, mv_Bidx);

mv_idx = sum(mm_GGNetMat,1)' > 0 & sum(mm_GGNetMat,2) > 0;
mc_Geneymbols = mc_ComGenes(mv_idx);

mn_D = length(mc_Geneymbols);

%- Match Genes in both sets 
mm_XGE = NaN*ones(mn_N, mn_D);
[mc_dummy4, mv_Aidx, mv_Bidx] = intersect(mc_Geneymbols, mc_gene_GE); 
mm_XGE(:,mv_Aidx) = mm_XGE_sub(:,mv_Bidx);

mm_GGINet = zeros(mn_D, mn_D);
[mc_dummy5, mv_Aidx, mv_Bidx] = intersect(mc_Geneymbols, mc_GGNet.gene); 
mm_GGINet(mv_Aidx,mv_Aidx) = mc_GGNet.Mat(mv_Bidx,mv_Bidx);

mm_GGINet = sparse(mm_GGINet);
% disp(['diagonal is zero:', find(diag(mm_GGINet)~=0)]);
% issymmetric(mm_GGINet)

%- Construct a pathway matrix (D X R)
if strcmpi(mstr_pathway,'pathdb')
    tmp = load('E:\BackUp\Code\NonnegativeMatrix\Mutation_05282013\data\gene_name_info_hg18');        
    hg18_gene_name = tmp.hg18_gene_name(:,1);
    
    tmp = load('E:\BackUp\Code\NonnegativeMatrix\Mutation_05282013\data\bipartite_PPI_module');
    module_idx_org = tmp.module_idx;
        
    % selection
    module_idx = module_idx_org;
    
    mc_PathwayInfo = cell(size(module_idx,1),1);
    for mn_i = 1:size(module_idx,1)
        mc_PathwayInfo{mn_i} = sprintf('%4dth pathway', mn_i);
    end
    
    mm_V0 = sparse([], [], [], mn_D, size(mc_PathwayInfo,1));%size(mc_PathwayInfo,1)    
    [mc_dummy6, mv_Aidx, mv_Bidx] = intersect(mc_Geneymbols, hg18_gene_name); 
    mm_V0(mv_Aidx,:) = module_idx(:,mv_Bidx)';
else
    if strcmpi(mstr_pathway, 'hallmark')
        tmp = load('./data/h.all.v6.2.symbols.gmt_member_Info.mat');
    else
        tmp = load(['./data/c2.cp.', mstr_pathway, '.v6.1.symbols.gmt_member_Info.mat']);
    end
    mc_PathwayInfo = tmp.mc_data.PathwayInfo;
    
    mm_V0 = sparse([], [], [], mn_D, size(mc_PathwayInfo,1));%size(mc_PathwayInfo,1)
    [mc_dummy6, mv_Aidx, mv_Bidx] = intersect(mc_Geneymbols, tmp.mc_data.genes); 
    mm_V0(mv_Aidx,:) = tmp.mc_data.G0(:,mv_Bidx)';
end

mn_R = size(mm_V0,2);

%--------------------------------------------------------------------------
% save('ROI_genes_Immnotherapy_repnose', 'mc_ROIgenes');
load('ROI_genes_Immnotherapy_repnose');

[mc_dummy, mv_Aidx, mv_Bidx] = intersect(mc_ROIgenes, mc_Geneymbols);

mm_XROI = NaN*ones(size(mm_XGE, 1), length(mc_ROIgenes));
mm_XROI(:, mv_Aidx) =  mm_XGE(:, mv_Bidx); 

mv_ResIDX = find(mm_U(:, 1) > 0);
mv_NonIDX = find(mm_U(:, 2) > 0);

mm_Results_GE = cell(length(mc_ROIgenes), 3); 
mm_Results_GE(:, 1) = mc_ROIgenes;
mm_Results_GE(:, 2) = cellstr(num2str(median(mm_XROI(mv_ResIDX, :), 1)'));
mm_Results_GE(:, 3) = cellstr(num2str(median(mm_XROI(mv_NonIDX, :), 1)'));

disp(median(mm_XROI(mv_ResIDX, :), 1))
disp(median(mm_XROI(mv_NonIDX, :), 1))
%--------------------------------------------------------------------------


%--- call BSNMTF
mstr_Linv_filename = ['MGA_', mstr_pathway];

mc_opt.MaxIters = 30;

mc_Ctypes = {'reponse', 'non-response'};
% out = myfunc_BSNMTF_ssp_final(mm_X, mm_U, mm_V0, mm_GGINet_mat, mc_opt);
mc_Solution = myfunc_BSNMTF_ssp_final...
                (mm_XGE, mm_U, mm_V0, mm_GGINet, mc_opt, mstr_Linv_filename);


save(['MGA_', mstr_pathway, '.mat'], 'mc_Solution', 'mc_PathwayInfo', 'mm_V0', 'mc_Geneymbols', '-v7.3');              


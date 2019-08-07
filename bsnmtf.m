function mc_Solution = bsnmtf(X, mm_U0, mm_Z0, mm_A, mc_opt)
    %----------------------------------------------------------------------
    % Bayesian semi-nonnegative matrix tri-factorization with 
    % spike-and-slab prior 
    %
    % (X = US\overline{V}^T = US(Z\cdotV)^T)
    % ---> p(V_jr, Z_jr_ G_{jr}) = N(V_jr*Z_jr|0,\sigma_{jr}^{V0})
    %                               \Phi(G_{jr})^{Z_{jr}}(1-\Phi(G_{jr}))^{1-Z_{jr}}
    % ---> g_r ~ GP(m_r, L), 
    %   where  m_r is constructed from the initial pathway membership info Z^0
    %          L is a laplace matrix constructed from PPI 
    % ---> variational inference
    % programed by Sunho Park (parks@ccf.org)
    %----------------------------------------------------------------------
    
    %--- setting
    [mn_N, mn_D] = size(X);
                
    mn_K = size(mm_U0,2);    
    mn_R = size(mm_Z0,2); % V \in R^{D X R}        
    
    [Iw, Jw, ~] = find(mm_Z0);
    V0_IDX = (Jw-1)*mn_D + Iw;
    
    % xi_+ and xi_- for the prior mean vectors of the GPs
    mr_mu_forzero = -5;
    mr_mu_forone = 5;
    
    mm_mprior = mr_mu_forzero*ones(mn_D, mn_R);
    mm_mprior(V0_IDX) = mr_mu_forone;
    
    % maximum iteration  
    mn_Maxiter = mc_opt.MaxIters;
    mflag_dispaly = 1;
        
    % for S
    mr_lambda = 10;
      
    % for L
    mn_L = myfunc_symNormalLap(mm_A);
    mm_Linv = mn_L\eye(mn_D,mn_D);
        
    mm_Linv_diag = full(diag(mm_Linv));
    
    % variational parameters
    global g_mu g_var Z_mu % mm_MatKRCHK
    
    global mm_MatDR01 mm_MatDR02 mm_MatDR03 mm_MatDR04 mm_MatDR05 mm_MatDR06
    global mm_MatND01 mv_VecD01
    global mm_MatKR01 mv_VecN01 mv_VecK01 mv_VecR01  
    
    %- X ~ US(ZoV'), where o is an element wise multiplication operator.
    % U
    mr_zeta = 1 - (1e-5);
    mr_1mzeta = 1 - mr_zeta;

    U(:,:) = mr_1mzeta + mr_zeta*mm_U0;
    U_db(:,:) = U.^2;
                
    mc_UIDX = cell(mn_K,1);
    for mn_k = 1:mn_K
        mc_UIDX{mn_k} = (mn_k-1)*mn_N+1:mn_k*mn_N;
    end
    
    % S
    S_mu = rand(mn_K, mn_R);     
    S_var = zeros(mn_K, mn_R);                       
    
    S_tmp_mu = zeros(mn_K, mn_R);               
    S_tmp_sig = zeros(mn_K, mn_R);                
    
    S_dbE = S_mu.^2;
 
    % V
    try
        [~, C] = kmeans(X, mn_R);
    catch
        mv_RandIDXSub = datasample(1:size(X,1), mn_R, 'Replace',true);
        
        C = X(mv_RandIDXSub,:);
    end
    
    V_mu  = C';
    V_var = 0.1*rand(mn_D, mn_R);   
    V_dbE = V_var + (V_mu.^2);
    
    V0_mu = 0;
    V0_sig = 1;
    V0_sig_inv = 1/V0_sig;            
    
    % [Zmu]_{jr} = p(Z_{jr}=1) = \rho_{jr}
    Z_mu = mm_Z0;    
      
    % \tau ~ Gamma(a0, b0)
    alpha_a0 = 0.1;
    alpha_a = alpha_a0 + 0.5*numel(X);
    alpha_b0 = 0.1;
        
    tau_mu = 0.1;
              
    % g_r ~ GP(0, L)
    g_mu = zeros(mn_D, mn_R);
    g_var = zeros(mn_D, mn_R);
    
    % setting for L-BFGS
    m_copt.maxIts = 20;
    m_copt.maxTotalIts = 200;
    m_copt.m = 5;
    m_copt.factor = 1e3;
    m_copt.pgtol = 1e-16;
    m = m_copt.m;
    
    mc_ParmIDX = cell(2,1);    
    mc_ParmIDX{1} = 1:(mn_R*mn_D);
    mc_ParmIDX{2} = ((mn_R*mn_D)+1):(2*mn_R*mn_D);
        
    l_par = -inf*ones(2*mn_R*mn_D,1);
    u_par =  inf*ones(2*mn_R*mn_D,1);
    nbd_par = isfinite(l_par) + isfinite(u_par) + 2*isinf(l_par).*isfinite(u_par);

    if ispc    
        nbd_par = int32(nbd_par);
    else        
        nbd_par = int64(nbd_par); 
    end
    
    factr   = m_copt.factor;
    pgtol   = m_copt.pgtol;
    maxIts  = m_copt.maxIts;
    maxTotalIts     = m_copt.maxTotalIts;
    printEvery = 0;
    iprint  = -1;  
    
   % matrics to hold intermediate values
    mm_MatKK01 = zeros(mn_K, mn_K); 
    mm_MatRR01 = zeros(mn_R, mn_R); 
       
    mm_MatDR01 = zeros(mn_D, mn_R); 
    mm_MatDR02 = zeros(mn_D, mn_R);
    mm_MatDR03 = zeros(mn_D, mn_R);
    mm_MatDR04 = zeros(mn_D, mn_R); 
    mm_MatDR05 = zeros(mn_D, mn_R);
    mm_MatDR06 = zeros(mn_D, mn_R);  
    
    mm_MatND01 = zeros(mn_N, mn_D);
        
    mv_VecD01 = zeros(mn_D, 1);
    
    mm_MatNR01 = zeros(mn_N, mn_R);
    mm_MatNR02 = zeros(mn_N, mn_R);
    mm_MatNR03 = zeros(mn_N, mn_R);
    mm_MatNR04 = zeros(mn_N, mn_R);
    
    mm_MatKR01 = zeros(mn_K, mn_R);
    mm_MatKR02 = zeros(mn_K, mn_R);
    mm_MatKR03 = zeros(mn_K, mn_R);
    
    mv_VecR01 = zeros(mn_R, 1);
    
    mv_VecN01 = zeros(mn_N, 1);
    mv_VecN02 = zeros(mn_N, 1);
    
    mv_VecK01 = zeros(mn_K, 1);
    
    mv_Krand = zeros(mn_K, 1);
    mv_Rrand = zeros(mn_R, 1);
    
    mv_IDX1K = 1:mn_K;
    mv_IDX1R = 1:mn_R;
    
    % initialize the variational parameters for G_jr
    mv_Zparams = 0.1*randn(2*mn_R*mn_D, 1);
    
    g_mu(:,:) = mm_mprior;
    mv_Zparams(mc_ParmIDX{1}) = reshape(g_mu, [mn_D*mn_R, 1]);

    g_var(:,:) = exp(reshape(mv_Zparams(mc_ParmIDX{2}), [mn_D, mn_R]));
    
    % initialize Z    
    Z_mu(:,:) = mm_Z0;
    
    %--- main
    fprintf('iter:: [Obj = ||X-UV(Z*V)||^2] \n' )
    
    for mn_iter = 1:mn_Maxiter                                        
        %- common
        mm_MatDR01(:,:) = Z_mu.*V_dbE;
        mm_MatDR02(:,:) = Z_mu.*V_mu;
        mm_MatDR03(:,:) = mm_MatDR02.^2;        

        %- update S   
        mv_VecR01(:) = sum(mm_MatDR01 - mm_MatDR03,1);        
        
        mm_MatKK01(:,:) = U'*U;
        mm_MatRR01(:,:) = mm_MatDR02'*mm_MatDR02;
        
        % sigdb^S
        S_dbE(:,:) = tau_mu*(sum(U_db,1)')*sum(mm_MatDR01,1);               
        S_tmp_sig(:,:) = sqrt(S_dbE);
        
        mm_MatKR01(:,:) = (U'*(X*mm_MatDR02));
        mm_MatKR02(:,:) = (mm_MatKK01*S_mu)*mm_MatRR01;
        mm_MatKR03(:,:) = sum(U_db,1)'*sum(mm_MatDR03,1);
        
        mm_MatNR01(:,:) = U*S_mu;
        
        mv_Krand(:) = randperm(mn_K);
        mv_Rrand(:) = randperm(mn_R);
        for mn_k_org = 1:mn_K
            mn_k = mv_Krand(mn_k_org);
            
            mv_VecN01(:) = U(:,mn_k);
            
            for mn_r_org = 1:mn_R 
                mn_r = mv_Rrand(mn_r_org);
                
                mr_sigdb = 1/S_dbE(mn_k,mn_r);
                
                mr_S = S_mu(mn_k, mn_r);
                mr_Ref01 = mm_MatKR02(mn_k,mn_r) - (mr_S*mm_MatKR03(mn_k,mn_r));
                
                mv_VecN02(:) = mm_MatNR01(:,mn_r) - mr_S*mv_VecN01;                
                
                mr_mu = - mr_lambda + tau_mu*( (mm_MatKR01(mn_k,mn_r) - mr_Ref01) ...
                                                - mv_VecR01(mn_r)*sum(mv_VecN01.*mv_VecN02) );                
                mr_mu = mr_sigdb*mr_mu;
                
                S_tmp_mu(mn_k,mn_r) = mr_mu; 
                
                % update the statistics (mean and variation)               
                mr_sig = sqrt(mr_sigdb);
                if (mr_mu < -30 * mr_sig)
                    % exp = 1./(abs(mu)*tau)
                    % var = (1./(abs(mu)*tau))**2
                    mr_mu_tilde = mr_sigdb/(abs(mr_mu)); 
                    mr_var_tilde = (mr_sigdb/abs(mr_mu))^2;
                else
                    mr_tmp01 = -(mr_mu/mr_sig);
                    mr_tmp02 = myfunc_lambda(mr_tmp01);

                    mr_mu_tilde  = mr_mu + mr_sig*mr_tmp02; 
                    mr_var_tilde = mr_sigdb*( 1 - (mr_tmp02*(mr_tmp02-mr_tmp01)) );            
                end
                
                if (mr_mu_tilde < 0) || isinf(mr_mu_tilde) || isnan(mr_mu_tilde)
                    mr_mu_tilde = 0;
                end
                
                if (mr_mu_tilde < 0) || isinf(mr_mu_tilde) || isnan(mr_mu_tilde)
                    mr_var_tilde = 0;
                end
                
                S_mu(mn_k,mn_r) = mr_mu_tilde;
                S_var(mn_k,mn_r) = mr_var_tilde;
                                
                % update M = U*S
                mm_MatNR01(:,mn_r) = mm_MatNR01(:,mn_r) + (mr_mu_tilde-mr_S)*mv_VecN01; 
                
                % update M = U'*(U*S*(Z*V)')*V                
                mm_MatKR02(:,:) = mm_MatKR02 ...
                    + (mr_mu_tilde-mr_S)*(mm_MatKK01(:,mn_k)*mm_MatRR01(mn_r,:));             
            end
            
            % update: to reduce numerical unstability
            mm_MatNR01(:,:) = U*S_mu;
            mm_MatKR02(:,:) = (mm_MatKK01*S_mu)*mm_MatRR01;
        end                                    
        S_dbE(:,:) = S_var + (S_mu.^2);        
        
        %- update V
        mm_MatNR01(:,:) = U*S_mu;
        mm_MatNR02(:,:) = mm_MatNR01.^2;
        
        V_var(:,:) = 1./( V0_sig_inv ...
                        + (tau_mu*repmat(sum(mm_MatNR02+(U_db*S_var),1),[mn_D,1])) );
              
        mm_MatDR01(:,:) = X'*mm_MatNR01;        
        mm_MatDR02(:,:) = Z_mu.*V_mu;
        
        mv_Rrand(:) = randperm(mn_R);
        for mn_r_org = 1:mn_R 
            mn_r = mv_Rrand(mn_r_org);
            
            mv_IDXIXsub = mv_IDX1R;
            mv_IDXIXsub(mn_r) = [];
            
            mv_VecR01(:,:) = mm_MatNR01'*mm_MatNR01(:,mn_r);
            
            V_mu(:,mn_r) = (V0_sig_inv*V0_mu) ...
                         + tau_mu*( mm_MatDR01(:,mn_r) ...
                                    - (mm_MatDR02(:,mv_IDXIXsub)*mv_VecR01(mv_IDXIXsub)) );                
             
            V_mu(:,mn_r) = V_var(:,mn_r).*V_mu(:,mn_r);
            
            % update Z_mu = \rho = p(Z=1)
            mv_VecD01(:) = g_mu(:,mn_r)./sqrt(1+g_var(:,mn_r));
             
            mv_VecD01(:) = logphi(mv_VecD01) - logphi(-mv_VecD01)...
                         + 0.5*((V_mu(:,mn_r).^2)./V_var(:,mn_r)) + 0.5*log(V_var(:,mn_r));
    
            Z_mu(:,mn_r) = 1./(1+exp(-mv_VecD01));
            
            % update M = Z.*V
            mm_MatDR02(:,mn_r) = Z_mu(:,mn_r).*V_mu(:,mn_r);
        end                
        V_dbE(:,:) = V_var + (V_mu.^2);
                
        %- tau   
        mm_MatDR01(:,:) = Z_mu.*V_dbE;
        mm_MatDR02(:,:) = Z_mu.*V_mu;
        mm_MatDR03(:,:) = mm_MatDR02.^2;
        
        % (sum_k sum_r)^2
        % mm_MatND01: A = X - USV'
        mm_MatND01(:,:) = X - ((U*S_mu)*mm_MatDR02');
        mr_sum1 = sum( mm_MatND01(:).^2 );
        
        % sum_k sum_r
        mm_MatNR01(:,:) = U_db*S_dbE;
        mm_MatNR02(:,:) = U_db*S_mu.^2;
        
        % sum_k sum_r sum_k'
        mm_MatNR03(:,:) = zeros(mn_N, mn_R);
        for mn_k = 1:mn_K
            mv_IDXIXsub = mv_IDX1K;
            mv_IDXIXsub(mn_k) = [];
            
            mm_MatNR04(:,:) = U(:,mv_IDXIXsub)*S_mu(mv_IDXIXsub,:);
            mm_MatNR04(:,:) = bsxfun(@times,mm_MatNR04,U(:,mn_k));
            mm_MatNR04(:,:) = bsxfun(@times,mm_MatNR04,S_mu(mn_k,:));
            
            mm_MatNR03(:,:) = mm_MatNR03 + mm_MatNR04;
        end        
        mm_MatNR01 = mm_MatNR01 + mm_MatNR03;
        mm_MatNR02 = mm_MatNR02 + mm_MatNR03;                        

        mr_sum2 = sum(sum(mm_MatNR01,1).*sum(mm_MatDR01,1));
        mr_sum3 = -sum(sum(mm_MatNR02,1).*sum(mm_MatDR03,1));

        alpha_b = alpha_b0 + 0.5*( mr_sum1 + mr_sum2 + mr_sum3 );
        
        tau_mu = alpha_a/alpha_b;

        % Optimize g_mu & g_sig      
        fcn_wrapper();
        
        callF = @(x)myfunc_ParamOpt(x, Z_mu, mm_Linv, mm_Linv_diag, mm_mprior, mc_ParmIDX);
        
        callF_wrapped = @(x, varargin) ...
                fcn_wrapper(callF, [], maxIts, printEvery, x, varargin{:});               
    
        x = mv_Zparams + 0;
        [~, mv_Zparams(:)] = lbfgsb_wrapper(m, x, l_par, u_par, nbd_par, ...
            callF_wrapped, factr, pgtol, iprint, maxIts, maxTotalIts);     
                
        g_mu(:,:) = reshape(mv_Zparams(mc_ParmIDX{1}), [mn_D, mn_R]);
        g_var(:,:) = exp(reshape(mv_Zparams(mc_ParmIDX{2}), [mn_D, mn_R]));                

        if mflag_dispaly > 0
            mm_MatDR01(:,:) = Z_mu.*V_dbE;
            mm_MatDR02(:,:) = Z_mu.*V_mu;
            mm_MatDR03(:,:) = mm_MatDR02.^2;

            %- (sum_k sum_r)^2
            % mm_MatND01: A = X - USV'
            mm_MatND01(:,:) = X - ((U*S_mu)*mm_MatDR02');
            mr_sum1 = sum( mm_MatND01(:).^2 );

            mr_diff_MF = mr_sum1/numel(X);
            
            fprintf('%3d:: ave. reconstruction error = %.3f \n', mn_iter, mr_diff_MF);
        end
    end
    
    mc_Solution.U = U;
    mc_Solution.S = S_mu;
    mc_Solution.V = V_mu;
    mc_Solution.Z = Z_mu;  
end

function [mr_fval, mv_grad] = myfunc_ParamOpt...
                (mv_theta, Z_mu, mm_Kinv, mv_Kinvdiag, mm_mprior, mc_ParmIDX)
      
    global g_mu g_var 
    global mm_MatDR01 mm_MatDR02 mm_MatDR03 mm_MatDR04 mm_MatDR05 mm_MatDR06 
    
    %- settings
    mn_LenParams = length(mv_theta);
    if nargout > 1
        mv_grad = zeros(mn_LenParams, 1);
    end
        
    [mn_D, mn_R] = size(Z_mu); 

    %- Params definitions
    g_mu(:,:) = reshape(mv_theta(mc_ParmIDX{1}), [mn_D, mn_R]);
    g_var(:,:) = reshape(exp(mv_theta(mc_ParmIDX{2})), [mn_D, mn_R]);
    
    %- int q(Z)q(G) log p(Z|G)p(G)
    mm_MatDR01(:,:) = 1 + g_var;
    mm_MatDR02(:,:) = g_mu./sqrt(mm_MatDR01);
    
    if nargout > 1
        [mm_MatDR03(:,:), mm_MatDR04(:,:)] = logphi(mm_MatDR02);
        [mm_MatDR05(:,:), mm_MatDR06(:,:)] = logphi(-mm_MatDR02);     
    else
        mm_MatDR03(:,:) = logphi(mm_MatDR02);
        mm_MatDR05(:,:) = logphi(-mm_MatDR02);     
    end
    
    mr_sumZjr = (Z_mu(:)'*mm_MatDR03(:)) + ((1-Z_mu(:))'*mm_MatDR05(:));
    
    % int q(G) log p(G)
    mm_MatDR02(:,:) = g_mu - mm_mprior;  
    mm_MatDR03(:,:) = (mm_Kinv*mm_MatDR02);
    
    mr_sum_hg = -0.5*( sum(mv_Kinvdiag'*g_var) + (mm_MatDR02(:)'*mm_MatDR03(:)) );
          
    % H(g)
    mr_sum_entro = 0.5*sum(mv_theta(mc_ParmIDX{2}));    
    
    % the output function values
    mr_fval = -(mr_sumZjr + mr_sum_hg + mr_sum_entro); 
            
    if nargout > 1
        %- mu: int q(Z)q(G) log p(Z|G)p(G)        
        mm_MatDR02(:,:) = ((Z_mu.*mm_MatDR04) - ((1-Z_mu).*mm_MatDR06))...
                        ./sqrt(mm_MatDR01);
        
        %- mu: int q(g) log p(g)
        mv_grad(mc_ParmIDX{1}) = reshape(mm_MatDR02 - mm_MatDR03,[mn_D*mn_R,1]);

        %- sigvar: int q(Z)q(G) log p(Z|G)p(G)        
        mm_MatDR02(:,:) = -(mm_MatDR02.*(g_mu./mm_MatDR01));
        
        %- sigvar: H(G) & int q(G) log p(G)
        mm_MatDR05(:,:) = mm_MatDR02 ...
                        + bsxfun(@minus,1./(g_var+eps), mv_Kinvdiag);
        
        mv_grad(mc_ParmIDX{2}) = reshape(mm_MatDR05,[mn_D*mn_R,1]);
        
        % final output
        mv_grad(mc_ParmIDX{2}) = mv_grad(mc_ParmIDX{2}).*reshape(g_var,[mn_D*mn_R,1]);

        mv_grad(:) = -mv_grad;
    end
end


function mm_val = myfunc_lambda(x)
    [~, mm_val] = logphi(-x);
end


function L = myfunc_symNormalLap(A)
    [m_nN, m_nD] = size(A);
    
    if m_nN ~= m_nD
       disp('The matrix is not square !'); 
       return;
    end
    
    [rows, cols, vals] = find(A); 
    
    mv_sum_diag = sqrt(sum(A,2));
    mv_Matsub = mv_sum_diag(rows).*mv_sum_diag(cols);
    
    L = speye(m_nD, m_nD) - sparse(rows, cols, vals./mv_Matsub, m_nD, m_nD);    
end


%--------------------------------------------------------------------------
% Safe computation of logphi(z) = log(normcdf(z)) and its derivatives
%                    dlogphi(z) = normpdf(x)/normcdf(x).
% The function is based on index 5725 in Hart et al. and gsl_sf_log_erfc_e.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2013-11-13.
function [lp,dlp,d2lp,d3lp] = logphi(z)
  z = real(z);                                 % support for real arguments only
  lp = zeros(size(z));                                         % allocate memory
  id1 = z.*z<0.0492;                                 % first case: close to zero
  lp0 = -z(id1)/sqrt(2*pi);
  c = [ 0.00048204; -0.00142906; 0.0013200243174; 0.0009461589032;
       -0.0045563339802; 0.00556964649138; 0.00125993961762116;
       -0.01621575378835404; 0.02629651521057465; -0.001829764677455021;
       2*(1-pi/3); (4-pi)/3; 1; 1];
  f = 0; for i=1:14, f = lp0.*(c(i)+f); end, lp(id1) = -2*f-log(2);
  id2 = z<-11.3137;                                    % second case: very small
  r = [ 1.2753666447299659525; 5.019049726784267463450;
        6.1602098531096305441; 7.409740605964741794425;
        2.9788656263939928886 ];
  q = [ 2.260528520767326969592;  9.3960340162350541504;
       12.048951927855129036034; 17.081440747466004316; 
        9.608965327192787870698;  3.3690752069827527677 ];
  num = 0.5641895835477550741; for i=1:5, num = -z(id2).*num/sqrt(2) + r(i); end
  den = 1.0;                   for i=1:6, den = -z(id2).*den/sqrt(2) + q(i); end
  e = num./den; lp(id2) = log(e/2) - z(id2).^2/2;
  id3 = ~id2 & ~id1; lp(id3) = log(erfc(-z(id3)/sqrt(2))/2);  % third case: rest
  if nargout>1                                        % compute first derivative
    dlp = zeros(size(z));                                      % allocate memory
    dlp( id2) = abs(den./num) * sqrt(2/pi); % strictly positive first derivative
    dlp(~id2) = exp(-z(~id2).*z(~id2)/2-lp(~id2))/sqrt(2*pi); % safe computation
    if nargout>2                                     % compute second derivative
      d2lp = -dlp.*abs(z+dlp);             % strictly negative second derivative
      if nargout>3                                    % compute third derivative
        d3lp = -d2lp.*abs(z+2*dlp)-dlp;     % strictly positive third derivative
      end
    end
  end
end

%--------------------------------------------------------------------------
function [f,g] = fcn_wrapper(callF, errFcn, maxIts, printEvery, x, varargin )
    persistent k history
    if isempty(k), k = 1; end
    if nargin==0
        % reset persistent variables and return information
%         if ~isempty(history) && ~isempty(k) 
%             printFcn(k,history);
%             f = history(1:k,:);
%         end
        history = [];
        k = [];
        return;
    end
    if isempty( history )
        width       = 0;
        if iscell( errFcn ), width = length(errFcn);
        elseif ~isempty(errFcn), width = 1; end
        width       = width + 2; % include fcn and norm(grad) as well
        history     = zeros( maxIts, width );
    end
    % Find function value and gradient:
    [f,g] = callF(x);
    if nargin > 5
        outerIter = varargin{1}+1;
        
        history(outerIter,1)    = f;
        history(outerIter,2)    = norm(g,Inf); % g is not projected
        if isa( errFcn, 'function_handle' )
            history(outerIter,3) = errFcn(x);
        elseif iscell( errFcn )
            for j = 1:length(errFcn)
                history(outer_count,j+2) = errFcn{j}(x);
            end
        end
        
        if outerIter > k
            % Display info from *previous* input
            % Since this may be called several times before outerIter
            % is actually updated
            %         fprintf('At iterate %5d, f(x)= %.2e, ||grad||_infty = %.2e [MATLAB]\n',...
            %             k,history(k,1),history(k,2) );
            if ~isinf(printEvery) && ~mod(k,printEvery)
                printFcn(k,history);
            end
            k = outerIter;
        end                
    end
end


function printFcn(k,history)
    fprintf('Iter %5d, f(x) = %2e, ||grad||_infty = %.2e', ...
        k, history(k,1), history(k,2) );
    for col = 3:size(history,2)
        fprintf(', %.2e', history(k,col) );
    end
    fprintf('\n');
end

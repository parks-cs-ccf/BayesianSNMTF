function BSNTMF_plot(md_X, md_U, md_Z0, md_A, md_X_hat, md_S_hat, md_Z_hat, md_V_hat, n_figure)
md_colormap = cool;
s_filetype = {'.png', '.eps', '.m'};
n_fontsize = 20;
n_angle = 45;

% X
figure();
imagesc(md_X);
title('$\bf{X}$', 'Interpreter', 'latex');
colormap(md_colormap), colorbar(); caxis([-5,5]);
ylabel('Samples'), xlabel('Genes')
set(gca, 'xtick', 0:100:800, 'xticklabelrotation', n_angle);
set(gca,'FontSize',n_fontsize);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Bayesian/X', s_filetype{1}]);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Bayesian/X', s_filetype{2}], 'epsc');
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Bayesian/X', s_filetype{3}]);

% U
figure();
imagesc(md_U);
title('$\bf{U}$', 'Interpreter', 'latex');
colormap(md_colormap);
xlabel('Subgroups'), ylabel('Samples');
set(gca, 'xtick', 1:4, 'xticklabel', {'A', 'B', 'C', 'D'});
set(gca,'FontSize',n_fontsize);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Bayesian/U', s_filetype{1}]);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Bayesian/U', s_filetype{2}], 'epsc');
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Bayesian/U', s_filetype{3}]);

% Z0
figure();
imagesc(md_Z0');
title('$\bf{Z^0}$', 'Interpreter', 'latex');
colormap(md_colormap);
ylabel('Pathways'), xlabel('Genes')
set(gca, 'ytick', 1:4);
set(gca, 'xtick', 0:100:800, 'xticklabelrotation', n_angle);
set(gca,'FontSize',n_fontsize);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Bayesian/Z^0', s_filetype{1}]);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Bayesian/Z^0', s_filetype{2}], 'epsc');
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Bayesian/Z^0', s_filetype{3}]);

% A
figure();
imagesc(md_A);
title('$\bf{A}$', 'Interpreter', 'latex');
colormap(md_colormap);
ylabel('Genes'), xlabel('Genes')
set(gca, 'xtick', 0:100:800, 'xticklabelrotation', n_angle);
set(gca, 'ytick', 0:100:800);
set(gca,'FontSize',n_fontsize);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Bayesian/A', s_filetype{1}]);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Bayesian/A', s_filetype{2}], 'epsc');
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Bayesian/A', s_filetype{3}]);

% X_hat
figure();
imagesc(md_X_hat);
title('$\bf{\hat{X}}$', 'Interpreter', 'latex');
colormap(md_colormap);
ylabel('Samples'), xlabel('Genes')
set(gca, 'xtick', 0:100:800, 'xticklabelrotation', n_angle);
set(gca,'FontSize',n_fontsize);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Bayesian/X_hat', s_filetype{1}]);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Bayesian/X_hat', s_filetype{2}], 'epsc');
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Bayesian/X_hat', s_filetype{3}]);

% S_hat
figure();
% heatmap(md_S_hat, {'1','2','3','4'}, {'A','B','C','D'}, '%0.2f', ...
%     'UseFigureColormap', true, 'Colorbar', false, 'UseLogColormap', true);
imagesc(md_S_hat);
title('$\bf{\hat{S}}$', 'Interpreter', 'latex');
colormap(md_colormap);
ylabel('Subgroups'), xlabel('Pathways')
set(gca, 'ytick', 1:4, 'yticklabel', {'A', 'B', 'C', 'D'});
set(gca,'FontSize',n_fontsize);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Bayesian/S_hat', s_filetype{1}]);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Bayesian/S_hat', s_filetype{2}], 'epsc');
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Bayesian/S_hat', s_filetype{3}]);

% Z_hat
figure();
imagesc(md_Z_hat');
title('$\bf{\hat{Z}}$', 'Interpreter', 'latex');
colormap(md_colormap);
ylabel('Pathways'), xlabel('Genes')
set(gca, 'ytick', 1:8);
set(gca, 'xtick', 0:100:800, 'xticklabelrotation', n_angle);
set(gca,'FontSize',n_fontsize);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Bayesian/Z_hat', s_filetype{1}]);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Bayesian/Z_hat', s_filetype{2}], 'epsc');
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Bayesian/Z_hat', s_filetype{3}]);

% V_hat
figure();
imagesc(md_V_hat');
title('$\bf{\hat{V}}$', 'Interpreter', 'latex');
colormap(md_colormap);
ylabel('Pathways'), xlabel('Genes');
set(gca, 'ytick', 1:8);
set(gca, 'xtick', 0:100:800, 'xticklabelrotation', n_angle);
set(gca,'FontSize',n_fontsize);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Bayesian/V_hat', s_filetype{1}]);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Bayesian/V_hat', s_filetype{2}], 'epsc');
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Bayesian/V_hat', s_filetype{3}]);

end
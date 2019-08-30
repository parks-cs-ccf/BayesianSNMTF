function NMTF_plot(md_X_new, md_U, md_V0_new, md_A_new, md_X_new_hat, md_S_hat, md_V_hat, n_figure)
md_colormap = cool;
s_filetype = {'.png', '.eps', '.m'};
n_fontsize = 20;
n_angle = 45;

% X_new
figure();
imagesc(md_X_new);
title('$\bf{X_{new}}$', 'Interpreter', 'latex');
colormap(md_colormap), colorbar(); caxis([0,5]);
ylabel('Samples'), xlabel('Genes')
set(gca, 'xtick', 0:200:1600, 'xticklabelrotation', n_angle);
set(gca,'FontSize',n_fontsize);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Folding/X', s_filetype{1}]);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Folding/X', s_filetype{2}], 'epsc');
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Folding/X', s_filetype{3}]);

% % U
figure();
imagesc(md_U);
title('$\bf{U}$', 'Interpreter', 'latex');
colormap(md_colormap);
xlabel('Subgroups'), ylabel('Samples');
set(gca, 'xtick', 1:4, 'xticklabel', {'A', 'B', 'C', 'D'});
set(gca,'FontSize',n_fontsize);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Folding/U', s_filetype{1}]);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Folding/U', s_filetype{2}], 'epsc');
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Folding/U', s_filetype{3}]);

% V0_new
figure();
imagesc(md_V0_new');
title('$\bf{V^0_{new}}$', 'Interpreter', 'latex');
colormap(md_colormap);
ylabel('Pathways'), xlabel('Genes')
set(gca, 'ytick', 1:4);
set(gca, 'xtick', 0:200:1600, 'xticklabelrotation', n_angle);
set(gca,'FontSize',n_fontsize);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Folding/Z^0', s_filetype{1}]);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Folding/Z^0', s_filetype{2}], 'epsc');
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Folding/Z^0', s_filetype{3}]);

% A_new
figure();
imagesc(md_A_new);
title('$\bf{A_{new}}$', 'Interpreter', 'latex');
colormap(md_colormap);
ylabel('Genes'), xlabel('Genes')
set(gca, 'xtick', 0:200:1600, 'xticklabelrotation', n_angle);
set(gca, 'ytick', 0:200:1600);
set(gca,'FontSize',n_fontsize);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Folding/A', s_filetype{1}]);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Folding/A', s_filetype{2}], 'epsc');
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Folding/A', s_filetype{3}]);

% X_new_hat
figure();
imagesc(md_X_new_hat);
title('$\bf{\hat{X}_{new}}$', 'Interpreter', 'latex');
colormap(md_colormap);
ylabel('Samples'), xlabel('Genes')
set(gca, 'xtick', 0:200:1600, 'xticklabelrotation', n_angle);
set(gca,'FontSize',n_fontsize);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Folding/X_hat', s_filetype{1}]);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Folding/X_hat', s_filetype{2}], 'epsc');
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Folding/X_hat', s_filetype{3}]);

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
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Folding/S_hat', s_filetype{1}]);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Folding/S_hat', s_filetype{2}], 'epsc');
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Folding/S_hat', s_filetype{3}]);

% V_new_hat
figure();
imagesc(md_V_hat');
title('$\bf{\hat{V}_{new}}$', 'Interpreter', 'latex');
colormap(md_colormap);
ylabel('Pathways'), xlabel('Genes');
set(gca, 'ytick', 1:8);
set(gca, 'xtick', 0:200:1600, 'xticklabelrotation', n_angle);
set(gca,'FontSize',n_fontsize);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Folding/V_hat', s_filetype{1}]);
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Folding/V_hat', s_filetype{2}], 'epsc');
saveas(gcf, ['./publication/supp_experiment_', num2str(n_figure), '/Folding/V_hat', s_filetype{3}]);

end
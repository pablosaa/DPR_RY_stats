% script to plot the results from the statistical analysis
% for DPRns, DPRans with classes Convective, Stratiform, and
% rainrate thresholds.
%

clear all;
close all;

SO(1) = load('DPR_RY_stats_RRth0.20.mat');
SO(2) = load('DPR_RY_stats_RRth0.50.mat');
%Wiqr = [1.5 5 Inf]';
N_test = length(SO(1).Wiqr);
Tiqr = {'97.0','99.3','100'};

%% ---------------------------------------------------------
%% For the Overall case:
% Plotting for correlation CorrRR, BIAS and RMSE
XLims = {[0.4 1], [-0.95 0], [0 4]};
figure(1)
set(gcf, 'PaperPositionMode', 'auto', 'Position', [548 909 1017 450]);
% Corr. coef.
ax(1) = subplot(1, 3, 1);
pcr{1}(:, 1) = plot(SO(1).CorrR, [1:N_test], '--^');
hold on;
pcr{1}(:, 2) = plot(SO(2).CorrR, [1:N_test], '-^');

%ldg = legend(pcr{1}(:,1), {'Thr=0.2 mm h^{-1}','Thr=0.5 mm h^{-1}'});
xlabel('Correlation coefficient');
ylabel('Data included [%]');

% BIAS:
ax(2) = subplot(1, 3, 2);
pcr{2}(:, 1) = plot(SO(1).BIAS, [1:N_test], '--^');
hold on;
pcr{2}(:, 2) = plot(SO(2).BIAS, [1:N_test], '-^');
xlabel({'(DPR - RY)','BIAS [mm h^{-1}]'});
ldg = legend(ax(2),{'DPR_{ns} @ RR\geq0.2  |  ', 'DPR_{ans} @ RR\geq0.2',...
										'DPR_{ns} @ RR\geq0.5  |  ','DPR_{ans} @ RR\geq0.5'},...
						 'Location','north','Orientation','horizontal','NumColumns',2);

% RMSE:
ax(3) = subplot(1, 3, 3);
pcr{3}(:, 1) = plot(SO(1).ubRMSE, [1:N_test], '--^');
hold on;
pcr{3}(:, 2) = plot(SO(2).ubRMSE, [1:N_test], '-^');
xlabel('ubRMSD [mm h^{-1}]');

%set(pcr, 'MarkerSize', 7, 'LineWidth', 2);
cellfun(@(x) set(x, 'MarkerSize', 7, 'LineWidth', 2, 'Color', [0.6 .7 0]), pcr);
cellfun(@(x) set(x(2,:), 'MarkerSize', 10, 'MarkerFaceColor', [1 1 1]*0.5,...
								 'Marker', 'o', 'Color', [.2 .5 .8]), pcr);

set(ax, 'YTick', [1:N_test], 'YLim', [0.5 3.5], 'XGrid', 'on', 'TickDir', 'out',...
		'PlotBoxAspectRatio', [0.9 0.8 0.8778], 'FontSize', 13);
set(ax(1),  'YTickLabel', Tiqr);
set(ax(2:3), 'YTickLabel', '');
tmp = get(ax, 'Position');
tmp = cellfun(@(x) [x(1) x(2)*1.05 x(3)*1.25 x(4)], tmp, 'UniformOutput', 0);
arrayfun(@(i) set(ax(i), 'Position', tmp{i}, 'XLim', XLims{i}), [1:3]);

%% -------------------------------------------------------------
%% For the precipitation type cases:

figure(2);
set(gcf, 'PaperPositionMode', 'auto');

iconv = [1 3];  % second colums is for DPR_ans
istra = [2 4];
% Corr. Coef.
bx(1,1) = subplot(3,2,1);
tyr{1}(:, 1) = plot( [1:N_test], SO(1).TypeCorrR(:, iconv), '--^');
hold on;
tyr{1}(:, 2) = plot( [1:N_test], SO(2).TypeCorrR(:, iconv), '-^');
ylabel('Correlation');
title('Convective');

bx(1,2) = subplot(3,2,2);
tyr{2}(:, 1) = plot( [1:N_test], SO(1).TypeCorrR(:, istra), '--^');
hold on;
tyr{2}(:, 2) = plot( [1:N_test], SO(2).TypeCorrR(:, istra), '-^');
title('Stratiform');

% BIAS:
bx(2,1) = subplot(3,2,3);
tyr{3}(:, 1) = plot([1:N_test], SO(1).TypeBIAS(:, iconv), '--^');
hold on;
tyr{3}(:, 2) = plot([1:N_test], SO(2).TypeBIAS(:, iconv), '-^');
ylabel('BIAS [mm hr^{-1}]');

bx(2,2) = subplot(3,2,4);
tyr{4}(:, 1) = plot([1:N_test], SO(1).TypeBIAS(:, istra), '--^');
hold on;
tyr{4}(:, 2) = plot([1:N_test], SO(2).TypeBIAS(:, istra), '-^');

% ubRMSE:
bx(3,1) = subplot(3,2,5);
tyr{5}(:, 1) = plot([1:N_test], SO(1).TypeubRMSE(:, iconv), '--^');
hold on;
tyr{5}(:, 2) = plot([1:N_test], SO(2).TypeubRMSE(:, iconv), '-^');
ylabel('ubRMSE [mm hr^{-1}]');
xlabel('Data included [%]');

bx(3,2) = subplot(3,2,6);
tyr{6}(:, 1) = plot([1:N_test], SO(1).TypeubRMSE(:, istra), '--^');
hold on;
tyr{6}(:, 2) = plot([1:N_test], SO(2).TypeubRMSE(:, istra), '-^');
xlabel('Data included [%]');

cellfun(@(x) set(x, 'MarkerSize', 7, 'LineWidth', 2, 'Color', [0.6 .7 0]), tyr);
cellfun(@(x) set(x(2,:), 'MarkerSize', 10, 'MarkerFaceColor', [1 1 1]*0.5,...
								 'Marker', 'o', 'Color', [.2 .5 .8]), tyr);

LimY(1,:) = [.5 .82];
LimY(2,:) = [-1.1 -0.3];
LimY(3,:) = [.8 2.4];

tmp0 = get(bx, 'Position');
tmp2 = cellfun(@(x) [1 0.97 1.18 1.27].*x, tmp0, 'UniformOutput', 0);

for i=1:3,
	set(bx(i,:), 'YLim', LimY(i,:));
	arrayfun(@(j) set(bx(i,j), 'Position', tmp2{i+3*(j-1)}), [1:2]);
end
set(bx(:,2), 'YTickLabel', '');
set(bx(1:2,:), 'XTickLabel', '');
set(bx(3,:), 'XTick', [1:3], 'XTickLabel', Tiqr);
set(bx, 'XLim', [.6 3.4], 'FontSize', 12, 'TickDir', 'out', 'YGrid', 'on');



% end of script.

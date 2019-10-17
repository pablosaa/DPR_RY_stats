% script for the statistics for Veli's paper on GPM and RADOLAN
%

clear all;
close all;

PLOT_FLAG = true;
MATF_FLAG = false;
PRINT_FLAG = false;

RRth = 0.5; % rain rate threshold to consider 0.2 (DPR specs), 0.5 mm/h (Technic)
Wiqr = [1.5 3 Inf]';  % sigma factors
N_test = 2; % number of test e.g. 2 for DPRns and DPRans.

load('../data/DPR_RADOLAN.mat');

% indices for DPRns and RY above threshold:
ii = ~(isnan(RY) | isnan(DPRns)) & (DPRns>=RRth & RY>=RRth) ;

% indices for DPRAns and RY above threshold:
kk = ~(isnan(RY) | isnan(DPRans)) & (DPRans>=RRth & RY>=RRth) & DPR_hip<1;

Xgr{1} = RY(ii);
Ysr{1} = double(DPRns(ii));

Xgr{2} = RY(kk);
Ysr{2} = double(DPRans(kk));

delRR = cellfun(@minus, Ysr, Xgr, 'UniformOutput', 0);
IQR   = cellfun(@(x) quantile(x, [.25 .75]), delRR, 'UniformOutput', 0);

% sigma times factor [N_wiqr x N_test]:
sigma = repmat(Wiqr, 1, N_test); %.*cellfun(@std, delRR);


for i=1:N_test,
	% Low and Up Threshold (column represent test):
	[MC Low_th{i} Ups_th{i}] = MedCouple(delRR{i}, Wiqr);
	%Low_th{i} = IQR{i}(:,1) - sigma(:,i).*diff(IQR{i});
	%Ups_th{i} = IQR{i}(:,2) + sigma(:,i).*diff(IQR{i});
	InIdx{i}  = arrayfun(@(a,b) find(delRR{i}>=a & delRR{i}<=b),...
											 Low_th{i}, Ups_th{i}, 'UniformOutput', 0);	

	% Correlation coefficients:
	TMPvar = cellfun(@(x) corrcoef(Xgr{i}(x), Ysr{i}(x)), InIdx{i},'UniformOutput',0);
	CorrR(:,i) = cellfun(@(x) x(1,2), TMPvar);

	% BIAS:
	BIAS(:,i) = cellfun(@(x) mean(Ysr{i}(x) - Xgr{i}(x)), InIdx{i});

	% RMSE:
	RMSE(:,i) = cellfun(@(x) sqrt(mean((Ysr{i}(x) - Xgr{i}(x)).^2)), InIdx{i});

	% Number of points:
	Ntot(:,i) = cellfun(@length, InIdx{i});

end

% ubRMSE:
ubRMSE = sqrt(RMSE.^2 - BIAS.^2);

% Percentage of data per class:
Nperc = Ntot./cellfun(@length, delRR);


%% Plotting the precipitation DPRns and DPRans database:
for i=1:2,
	if ~PLOT_FLAG, continue; end
	idxstat = 1;
	RRmax = 1e3;
	NNmax = 1e3;
	figure(i);
	set(gcf, 'PaperPositionMode','auto', 'Position',[611 508 678 507]);

	% 1:1 line:
	RRx = logspace(log10(RRth), log10(RRmax), 100);
	
	% Limits line:
	delRRmean = median(delRR{i});
	RRy = RRx + delRRmean;
	RRlow = Low_th{i}(idxstat,:)+RRy;
	RRhig = Ups_th{i}(idxstat,:)+RRy;

	% best fit line:
	pfit = polyfit(log10(Xgr{i}(InIdx{i}{idxstat})), log10(Ysr{i}(InIdx{i}{idxstat})),1);
	Yfit = polyval(pfit, RRx);

	% 2D histogram:
	XEDGES = logspace(log10(RRth), log10(RRmax), 70);
	NN = histcounts2(Xgr{i}, Ysr{i}, XEDGES, XEDGES);
	NN(NN<1) = NaN;
	pcolor(XEDGES(2:end),XEDGES(2:end),log10(NN'));
	shading flat;
	xlabel('RADOLAN RY [mm h^{-1}]');
	if i==1,
		ylabel('DPRns [mm h^{-1}]');
	else
		ylabel('DPRans [mm h^{-1}]');
	end
	txth = text(1.2*RRth, RRmax/2, sprintf('Rainrate Threshold\nRR_{thr} \\geq %3.1f [mm h^{-1}]', RRth),...
						 'FontSize',14);
	hbar = colorbar;
	ylabel(hbar, 'log10(N)');
	colormap((summer + jet)/2);
	hold on;
	t11 = plot(RRx, RRx, '--k', 'LineWidth', 2);

	l99 = loglog([RRx 1 RRx],[RRlow NaN RRhig],'--r','LineWidth',2);
	lfit = loglog(RRx, Yfit, '-', 'LineWidth', 2, 'Color', [0 .447 .8710]);
	set(gca,'XScale','log','YScale','log','TickDir','out',...
			'XLim',[RRth RRmax],'YLim',[RRth RRmax], 'TickLength',[0.02 0.04],...
			'CLim', [0 log10(NNmax)], 'FontSize',13);
	set(hbar, 'FontSize',11, 'TickDir', 'out');
	% legend
	txl = text(50, RRmax/2, sprintf('R = %4.3f\nBIAS = %4.2f\nRMSD = %3.2f\nubRMSD = %3.2f\nN = %d',...
															CorrR(idxstat,i), BIAS(idxstat,i), RMSE(idxstat,i),...
															ubRMSE(idxstat,i), Ntot(idxstat,i)),...
						 'FontSize',13, 'BackgroundColor', 'w', 'EdgeColor', [.4 .4 .4]);
	lgd = legend([t11 l99 lfit], {'1:1', sprintf('%3.1f%%',100*Nperc(idxstat,i)), 'FIT'},'Location','southeast');
end


%% ----------------------------------------------------------------
%% ----------------------------------------------------------------
%% FOR Precipitation Type (convective, stratiform)
N_type = 4; 
type_flag{1} = round(DPRns_ty/1e7)==1 & ii;  % 1=convective DPRns
type_flag{2} = round(DPRns_ty/1e7)==2 & ii;  % 2=stratiform DPRns
type_flag{3} = round(DPRns_ty/1e7)==1 & kk;  % 1=convective DPRans
type_flag{4} = round(DPRns_ty/1e7)==2 & kk;  % 2=stratiform DPRans


TypXgr{1} = RY(type_flag{1});
TypYsr{1} = double(DPRns(type_flag{1}));

TypXgr{2} = RY(type_flag{2});
TypYsr{2} = double(DPRns(type_flag{2}));

TypXgr{3} = RY(type_flag{3});
TypYsr{3} = double(DPRans(type_flag{3}));

TypXgr{4} = RY(type_flag{4});
TypYsr{4} = double(DPRans(type_flag{4}));


TypdelRR = cellfun(@minus, TypYsr, TypXgr, 'UniformOutput', 0);
IQR   = cellfun(@(x) quantile(x,[.25 .75]), TypdelRR, 'UniformOutput', 0);

% sigma times factor [N_wiqr x N_test]:
sigma = repmat(Wiqr, 1, N_type); %.*cellfun(@std, TypdelRR);

for i=1:N_type,
	% Low and Up Threshold (column represent type):
	[MC Low_th{i} Ups_th{i}] = MedCouple(TypdelRR{i}, Wiqr);
	%Low_th{i} = IQR{i}(:,1) - sigma(:,i).*diff(IQR{i});
	%Ups_th{i} = IQR{i}(:,2) + sigma(:,i).*diff(IQR{i});
	
	InIdx{i}  = arrayfun(@(a,b) find(TypdelRR{i}>=a & TypdelRR{i}<=b),...
											 Low_th{i}, Ups_th{i}, 'UniformOutput', 0);	

	
	% Correlation coefficients:
	TMPvar = cellfun(@(x) corrcoef(TypXgr{i}(x), TypYsr{i}(x)), InIdx{i},'UniformOutput',0);
	TypeCorrR(:,i) = cellfun(@(x) x(1,2), TMPvar);

	% BIAS:
	TypeBIAS(:,i) = cellfun(@(x) mean(TypYsr{i}(x) - TypXgr{i}(x)), InIdx{i});

	% RMSE:
	TypeRMSE(:,i) = cellfun(@(x) sqrt(mean((TypYsr{i}(x) - TypXgr{i}(x)).^2)), InIdx{i});

	% Number of points:
	TypeNtot(:,i) = cellfun(@length, InIdx{i});
end
% ubRMSE:
TypeubRMSE = sqrt(TypeRMSE.^2 - TypeBIAS.^2);

% Percentage of data per class:
TypeNperc = TypeNtot./cellfun(@length, TypdelRR);

%% Plotting the precipitation type database:
for i=1:N_type,
	if ~PLOT_FLAG, continue; end
	NNmax = [2 2.5 2 2.5]; %0.5e3;
	figure(i+2);
	delRRmean = mean(TypdelRR{i});
	RRx = logspace(log10(RRth), log10(RRmax), 100);
	RRy = RRx + delRRmean;

	% best fit line:
	Typfit = polyfit(log10(TypXgr{i}(InIdx{i}{idxstat})), log10(TypYsr{i}(InIdx{i}{idxstat})),1);
	TyYfit = polyval(Typfit, RRx);

	XEDGES = logspace(log10(RRth), 2, 70);
	NN = histcounts2(TypXgr{i}, TypYsr{i}, XEDGES, XEDGES);
	NN(NN<1) = NaN;
	pcolor(XEDGES(2:end), XEDGES(2:end), log10(NN'));
	shading flat;
	colorbar
	colormap((summer + jet)/2);
	xlabel('RADOLAN RY [mm h^{-1}]');
	if i<3,
		ylabel('DPRns [mm h^{-1}]');
	else
		ylabel('DPRans [mm h^{-1}]');
	end

	hold on;
	t11 = plot(RRx, RRx, '--k','LineWidth',2);
	%m50 = plot(RRx, RRy, '--b');
	l99 = loglog([RRx 1 RRx],[Low_th{i}(2,:)+RRy NaN Ups_th{i}(2,:)+RRy],'--r','LineWidth',2);
	lfit = loglog(RRx, TyYfit, '-', 'LineWidth', 2, 'Color', [0 .447 .8710]);
	set(gca,'XScale','log','YScale','log','TickDir','out', 'CLim', [0 NNmax(i)],...
			'XLim',[RRth RRmax*(1-RRth)],'YLim',[RRth RRmax*(1-RRth)],'TickLength',[0.02 0.04]);

	% legend
	txl = text(50, 2, sprintf('R = %4.3f\nBIAS = %4.2f\nRMSD = %3.2f\nubRMSD = %3.2f\nN = %d',...
														TypeCorrR(idxstat,i), TypeBIAS(idxstat,i), TypeRMSE(idxstat,i),...
														TypeubRMSE(idxstat,i), TypeNtot(idxstat,i)),...
						 'FontSize',13, 'BackgroundColor', 'w', 'EdgeColor', [.4 .4 .4]);
	lgd = legend([t11 l99 lfit], {'1:1', sprintf('%3.1f%%',100*TypeNperc(idxstat,i)), 'FIT'},'Location','northwest');

end

%% ---------------------------------
%% Storing the results in a mat file
%%

if MATF_FLAG,
	save(sprintf('../data/DPR_RY_stats_RRth%03.2f.mat',RRth),'CorrR','BIAS','RMSE','ubRMSE','Ntot','Nperc',...
			 'TypeCorrR','TypeBIAS','TypeRMSE','TypeubRMSE','TypeNtot','TypeNperc','RRth','Wiqr');
end

if PRINT_FLAG,
	print('-f1', '-dpng', sprintf('./plots/DPRns_RY_th%02.1f.png', RRth) );
	print('-f2', '-dpng', sprintf('./plots/DPRans_RY_th%02.1f.png', RRth) );
	print('-f3', '-dpng', sprintf('./plots/Convective_DPRns_RY_th%02.1f.png', RRth) );
	print('-f4', '-dpng', sprintf('../plots/Stratiform_DPRns_RY_th%02.1f.png', RRth) );
	print('-f5', '-dpng', sprintf('../plots/Convective_DPRans_RY_th%02.1f.png', RRth) );
	print('-f6', '-dpng', sprintf('../plots/Stratiform_DPRans_RY_th%02.1f.png', RRth) );

end

return;

%% TO BE FINISHED WHEN THE DATA WITH FLAG_PHASE for DPR_ans COME:
%% --------------------------------------------------------------------
%% FOR Hydrometeor Phase (liquid, solid, mixed)  !!! NOT VALID, new database to come
N_phase = 2;
phase_flag{1} = round(DPRns_ph/1e2)==1 & kk; % 1=mixed?
phase_flag{2} = round(DPRns_ph/1e2)==2 & kk; % 2=liquid

Xgr{1} = RY(phase_flag{1});
Ysr{1} = double(DPRans(phase_flag{1}));

Xgr{2} = RY(phase_flag{2});
Ysr{2} = double(DPRans(phase_flag{2}));

delRR = cellfun(@minus, Ysr, Xgr, 'UniformOutput', 0);
IQR   = cellfun(@(x) quantile(x,[.25 .75]), delRR, 'UniformOutput', 0);

% sigma times factor [N_wiqr x N_test]:
sigma = Wiqr.*cellfun(@std, delRR);

for i=1:N_phase,
	% Low and Up Threshold (column represent type):
	Low_th{i} = IQR{i}(:,1) - sigma(:,i).*diff(IQR{i});
	Ups_th{i} = IQR{i}(:,2) + sigma(:,i).*diff(IQR{i});
	InIdx{i}  = arrayfun(@(a,b) find(delRR{i}>=a & delRR{i}<=b),...
									 Low_th{i}, Ups_th{i}, 'UniformOutput', 0);	

	% Correlation coefficients:
	TMPvar = cellfun(@(x) corrcoef(Xgr{i}(x), Ysr{i}(x)), InIdx{i},'UniformOutput',0);
	PhaseCorrR(:,i) = cellfun(@(x) x(1,2), TMPvar);

	% BIAS:
	PhaseBIAS(:,i) = cellfun(@(x) mean(Ysr{i}(x) - Xgr{i}(x)), InIdx{i});

	% RMSE:
	PhaseRMSE(:,i) = cellfun(@(x) sqrt(mean((Ysr{i}(x) - Xgr{i}(x)).^2)), InIdx{i});

	% Total Number of points per class:
	PhaseNtot(:,i) = cellfun(@length, InIdx{i});
end
% ubRMSE:
PhaseubRMSE = sqrt(PhaseRMSE.^2 - PhaseBIAS.^2);

% end of script

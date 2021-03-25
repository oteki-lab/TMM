clearvars
close all

%% Import the data
[~, ~, raw] = xlsread('PL210312.xlsx','III201116A','A2:K513');

%% Create output variable
data = reshape([raw{:}],size(raw));

%% Allocate imported array to column variable names
wavelength = data(:,1);
% wavelength = data(:,1);
ref_10 = data(:,2);
ref_20 = data(:,3);
ref_30 = data(:,4);
ref_40 = data(:,5);
ref_50 = data(:,6);
ref_60 = data(:,7);
ref_70 = data(:,8);
ref_80 = data(:,9);
ref_90 = data(:,10);
ref_100 = data(:,11);

clear raw data

%% Import the data
[~, ~, raw] = xlsread('PL210312.xlsx','III201116A1','A2:K513');

%% Create output variable
data = reshape([raw{:}],size(raw));

%% Allocate imported array to column variable names
% wavelength = data(:,1);
FP_10 = data(:,2);
FP_20 = data(:,3);
FP_30 = data(:,4);
FP_40 = data(:,5);
FP_50 = data(:,6);
FP_60 = data(:,7);
FP_70 = data(:,8);
FP_80 = data(:,9);
FP_90 = data(:,10);
FP_100 = data(:,11);

clear raw data

%% PL ratios

enhancement10=FP_10./ref_10;
enhancement20=FP_20./ref_20;
enhancement30=FP_30./ref_30;
enhancement40=FP_40./ref_40;
enhancement50=FP_50./ref_50;
enhancement60=FP_60./ref_60;
enhancement70=FP_70./ref_70;
enhancement80=FP_80./ref_80;
enhancement90=FP_90./ref_90;
enhancement100=FP_100./ref_100;

PL_ratio_ref_20=ref_20./ref_20;
PL_ratio_ref_40=ref_40./ref_20;
PL_ratio_ref_60=ref_60./ref_20;
PL_ratio_ref_80=ref_80./ref_20;
PL_ratio_ref_100=ref_100./ref_20;

PL_ratio_FP_20=FP_20./FP_20;
PL_ratio_FP_40=FP_40./FP_20;
PL_ratio_FP_60=FP_60./FP_20;
PL_ratio_FP_80=FP_80./FP_20;
PL_ratio_FP_100=FP_100./FP_20;

%% Compare with theory

load('Absorption_enhancement_theory')
lambda=lambda*1e9;
lambda_min=min(lambda);
lambda_max=max(lambda);

%% Plots

%PL spectra
figure
plot(wavelength,ref_50,'linewidth',3)
hold on
plot(wavelength,FP_50,'linewidth',3)
xlabel('$\lambda \: (\mathrm{nm})$', 'Interpreter', 'LaTeX')
ylabel('Photoluminescence (arb. units)','Interpreter','Latex')
set(gca,'Fontsize',16)
set(gca,'XMinorTick','on','YMinorTick','on','YScale','log')
set(gcf,'color','w');
box on
legend({'before', 'after'},'Interpreter','Latex')
xlim([lambda_min, lambda_max])
ylim([1 1e5])
yticks([1 1e1 1e2 1e3 1e4 1e5])

% %PL spectra
% figure
% plot(wavelength,ref_20,'linewidth',3)
% hold on
% plot(wavelength,ref_40,'linewidth',3)
% plot(wavelength,ref_60,'linewidth',3)
% plot(wavelength,ref_80,'linewidth',3)
% plot(wavelength,ref_100,'linewidth',3)
% xlabel('Wavelength (nm)', 'Interpreter', 'LaTeX')
% ylabel('Photoluminescence (arb. units)','Interpreter','Latex')
% set(gca,'Fontsize',16)
% set(gca,'XMinorTick','on','YMinorTick','on','YScale','log')
% set(gcf,'color','w');
% box on
% legend({'20\% power', '40\% power', '60\% power', '80\% power', '100\% power'},'Interpreter','Latex')
% xlim([lambda_min, lambda_max])
% ylim([1 1e5])
% yticks([1 1e1 1e2 1e3 1e4 1e5])
% 
% figure
% plot(wavelength,FP_20,'linewidth',3)
% hold on
% plot(wavelength,FP_40,'linewidth',3)
% plot(wavelength,FP_60,'linewidth',3)
% plot(wavelength,FP_80,'linewidth',3)
% plot(wavelength,FP_100,'linewidth',3)
% xlabel('Wavelength (nm)', 'Interpreter', 'LaTeX')
% ylabel('Photoluminescence (arb. units)','Interpreter','Latex')
% set(gca,'Fontsize',16)
% set(gca,'XMinorTick','on','YMinorTick','on','YScale','log')
% set(gcf,'color','w');
% box on
% legend({'20\% power', '40\% power', '60\% power', '80\% power', '100\% power'},'Interpreter','Latex')
% xlim([lambda_min, lambda_max])
% ylim([1 1e5])
% yticks([1 1e1 1e2 1e3 1e4 1e5])
% 
% 
% 
% 
% %PL ratios
% figure
% plot(wavelength,PL_ratio_ref_20,'linewidth',2)
% hold on
% plot(wavelength,PL_ratio_ref_40,'linewidth',2)
% plot(wavelength,PL_ratio_ref_60,'linewidth',2)
% plot(wavelength,PL_ratio_ref_80,'linewidth',2)
% plot(wavelength,PL_ratio_ref_100,'linewidth',2)
% xlabel('Wavelength (nm)', 'Interpreter', 'LaTeX')
% ylabel('Photoluminescence (arb. units)','Interpreter','Latex')
% set(gca,'Fontsize',16)
% set(gca,'XMinorTick','on','YMinorTick','on')
% set(gcf,'color','w');
% box on
% legend({'20\% power', '40\% power', '60\% power', '80\% power', '100\% power'},'Interpreter','Latex')
% xlim([lambda_min, lambda_max])
% 
% figure
% plot(wavelength,PL_ratio_FP_20,'linewidth',2)
% hold on
% plot(wavelength,PL_ratio_FP_40,'linewidth',2)
% plot(wavelength,PL_ratio_FP_60,'linewidth',2)
% plot(wavelength,PL_ratio_FP_80,'linewidth',2)
% plot(wavelength,PL_ratio_FP_100,'linewidth',2)
% xlabel('Wavelength (nm)', 'Interpreter', 'LaTeX')
% ylabel('Photoluminescence (arb. units)','Interpreter','Latex')
% set(gca,'Fontsize',16)
% set(gca,'XMinorTick','on','YMinorTick','on')
% set(gcf,'color','w');
% box on
% legend({'20\% power', '40\% power', '60\% power', '80\% power', '100\% power'},'Interpreter','Latex')
% xlim([lambda_min, lambda_max])
% 

%Absorption enhancement
figure
plot(wavelength,enhancement10,'linewidth',2)
hold on
plot(wavelength,enhancement20,'linewidth',2)
plot(wavelength,enhancement30,'linewidth',2)
plot(wavelength,enhancement40,'linewidth',2)
plot(wavelength,enhancement50,'linewidth',2)
plot(wavelength,enhancement60,'linewidth',2)
plot(wavelength,enhancement70,'linewidth',2)
plot(wavelength,enhancement80,'linewidth',2)
plot(wavelength,enhancement90,'linewidth',2)
plot(wavelength,enhancement100,'linewidth',2)
xlabel('Wavelength (nm)', 'Interpreter', 'LaTeX')
ylabel('Photoluminescence ratio','Interpreter','Latex')
set(gca,'Fontsize',16)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w');
box on
%legend({'20\% power', '40\% power', '60\% power', '80\% power', '100\% power'},'Interpreter','Latex')
xlim([lambda_min, lambda_max])
ylim([0 5])


figure
%plot(wavelength,enhancement20/max(enhancement20(50:200)),'linewidth',2)
yyaxis left
plot(lambda,I_mean,'--','linewidth',2,'color',[0.00 0.45 0.74])
ylabel('$F(\lambda)$','Interpreter','Latex')
ylim([0 15])
%hold on
%plot(wavelength,enhancement50/max(enhancement50(50:200)),'linewidth',2)
%plot(wavelength,enhancement60/max(enhancement60(50:200)),'linewidth',2)
%plot(wavelength,enhancement80/max(enhancement80(50:200)),'linewidth',2)
%plot(wavelength,enhancement100/max(enhancement100(50:200)),'linewidth',2)
yyaxis right
plot(wavelength,enhancement50,'linewidth',2,'color',[1 0 0])
xlabel('$\lambda \: (\mathrm{nm})$', 'Interpreter', 'LaTeX')
ylabel('PL ratio','Interpreter','Latex')
set(gca,'Fontsize',16)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w');
ylim([0 3])
box on
xlim([lambda_min, lambda_max])
%legend({'Experiment20','Experiment40','Experiment50','Experiment60','Experiment80','Experiment100','Theory'})
%legend({'Experiment50','Theory'})

% 
% figure
% plot(wavelength,10*test,'linewidth',2)
% hold on
% plot(lambda,I_mean,'--','linewidth',2)
% xlabel('Wavelength (nm)', 'Interpreter', 'LaTeX')
% ylabel('Mean absorption enhancement','Interpreter','Latex')
% set(gca,'Fontsize',16)
% set(gca,'XMinorTick','on','YMinorTick','on')
% set(gcf,'color','w');
% box on
% xlim([lambda_min, lambda_max])
% legend({'Experiment','Theory'})
clearvars
close all

%% Import the data
[~, ~, raw] = xlsread('PL210107.xlsx','III201116A','A2:F513');

%% Create output variable
data = reshape([raw{:}],size(raw));

%% Allocate imported array to column variable names
wavelength = data(:,1);
ref_20 = data(:,2);
ref_40 = data(:,3);
ref_60 = data(:,4);
ref_80 = data(:,5);
ref_100 = data(:,6);

clear raw data

%% Import the data
[~, ~, raw] = xlsread('PL210107.xlsx','III201116A1','A2:K513');

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

enhancement20=FP_20./ref_20;
enhancement40=FP_40./ref_40;
enhancement60=FP_60./ref_60;
enhancement80=FP_80./ref_80;
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
plot(wavelength,ref_20,'linewidth',3)
hold on
plot(wavelength,ref_40,'linewidth',3)
plot(wavelength,ref_60,'linewidth',3)
plot(wavelength,ref_80,'linewidth',3)
plot(wavelength,ref_100,'linewidth',3)
xlabel('Wavelength (nm)', 'Interpreter', 'LaTeX')
ylabel('Photoluminescence (arb. units)','Interpreter','Latex')
set(gca,'Fontsize',16)
set(gca,'XMinorTick','on','YMinorTick','on','YScale','log')
set(gcf,'color','w');
box on
legend({'20\% power', '40\% power', '60\% power', '80\% power', '100\% power'},'Interpreter','Latex')
xlim([lambda_min, lambda_max])
ylim([1 1e5])
yticks([1 1e1 1e2 1e3 1e4 1e5])

figure
plot(wavelength,FP_20,'linewidth',3)
hold on
plot(wavelength,FP_40,'linewidth',3)
plot(wavelength,FP_60,'linewidth',3)
plot(wavelength,FP_80,'linewidth',3)
plot(wavelength,FP_100,'linewidth',3)
xlabel('Wavelength (nm)', 'Interpreter', 'LaTeX')
ylabel('Photoluminescence (arb. units)','Interpreter','Latex')
set(gca,'Fontsize',16)
set(gca,'XMinorTick','on','YMinorTick','on','YScale','log')
set(gcf,'color','w');
box on
legend({'20\% power', '40\% power', '60\% power', '80\% power', '100\% power'},'Interpreter','Latex')
xlim([lambda_min, lambda_max])
ylim([1 1e5])
yticks([1 1e1 1e2 1e3 1e4 1e5])




%PL ratios
figure
plot(wavelength,PL_ratio_ref_20,'linewidth',2)
hold on
plot(wavelength,PL_ratio_ref_40,'linewidth',2)
plot(wavelength,PL_ratio_ref_60,'linewidth',2)
plot(wavelength,PL_ratio_ref_80,'linewidth',2)
plot(wavelength,PL_ratio_ref_100,'linewidth',2)
xlabel('Wavelength (nm)', 'Interpreter', 'LaTeX')
ylabel('Photoluminescence (arb. units)','Interpreter','Latex')
set(gca,'Fontsize',16)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w');
box on
legend({'20\% power', '40\% power', '60\% power', '80\% power', '100\% power'},'Interpreter','Latex')
xlim([lambda_min, lambda_max])

figure
plot(wavelength,PL_ratio_FP_20,'linewidth',2)
hold on
plot(wavelength,PL_ratio_FP_40,'linewidth',2)
plot(wavelength,PL_ratio_FP_60,'linewidth',2)
plot(wavelength,PL_ratio_FP_80,'linewidth',2)
plot(wavelength,PL_ratio_FP_100,'linewidth',2)
xlabel('Wavelength (nm)', 'Interpreter', 'LaTeX')
ylabel('Photoluminescence (arb. units)','Interpreter','Latex')
set(gca,'Fontsize',16)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w');
box on
legend({'20\% power', '40\% power', '60\% power', '80\% power', '100\% power'},'Interpreter','Latex')
xlim([lambda_min, lambda_max])


%Absorption enhancement
figure
plot(wavelength,enhancement20,'linewidth',2)
hold on
plot(wavelength,enhancement40,'linewidth',2)
plot(wavelength,enhancement60,'linewidth',2)
plot(wavelength,enhancement80,'linewidth',2)
plot(wavelength,enhancement100,'linewidth',2)
xlabel('Wavelength (nm)', 'Interpreter', 'LaTeX')
ylabel('Photoluminescence ratio','Interpreter','Latex')
set(gca,'Fontsize',16)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w');
box on
legend({'20\% power', '40\% power', '60\% power', '80\% power', '100\% power'},'Interpreter','Latex')
xlim([lambda_min, lambda_max])


figure
plot(wavelength,10*enhancement20,'linewidth',2)
hold on
plot(wavelength,12*enhancement40,'linewidth',2)
plot(wavelength,13*enhancement60,'linewidth',2)
plot(wavelength,14*enhancement80,'linewidth',2)
plot(wavelength,15*enhancement100,'linewidth',2)
plot(lambda,I_mean,'--','linewidth',2)
xlabel('Wavelength (nm)', 'Interpreter', 'LaTeX')
ylabel('Mean absorption enhancement','Interpreter','Latex')
set(gca,'Fontsize',16)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w');
box on
xlim([lambda_min, lambda_max])
legend({'Experiment20','Experiment40','Experiment60','Experiment80','Experiment100','Theory'})

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
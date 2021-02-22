clear all

lambda=(0.920:0.001:1.340)*1e-6; %Wavelength vector
nb_lambda=width(lambda); %Number of wavelength

cal_abs=1; % calculate the absorption and field intensity in each layer
cal_field=1; % calculate the field intensity in each layer
transform_E=0; % Transform the data into a function of energy

back_b = 50;
surf_b = 25;
spacer = 25;
intensity = 12;

%% List of layers

n_0=1*ones(1,nb_lambda); %Incident medium
n_end=retindice_metal(lambda*1e6,1.9); %Last medium

% list of indices of the stack
n(1,:)=retindice_semicond(lambda*1e6,40);

% list of thicknesses of the stack
d(1)=1800;

[I, ~, stack]=TMM(nb_lambda,lambda,cal_abs,cal_field,transform_E,n_0,n_end,n,d);

I_mean = I(:,abs(lambda*1e6-1.192) < 0.001);    %I_mean = mean(I,2); % mean I for each hight
I_table = cat(2, I_mean, stack');
I_sorted = sortrows(I_table(back_b:end-surf_b,:),'descend');

QD_table = I_sorted(1,:);
for i=2:height(I_sorted)
    if abs(I_sorted(i,2)-QD_table(:,2))>spacer & I_sorted(i,1)>=intensity
        QD_table = cat(1, QD_table, I_sorted(i,:));
    end
end
QD_table = sortrows(QD_table,2);

%Plot
figure
plot(I_table(:,1),I_table(:,2),'Linewidth',3);
hold all
for i=1:height(QD_table)
    plot([0 QD_table(i,1)],[QD_table(i,2) QD_table(i,2)],'linewidth',1,'color','r')
end
ylim([0 sum(d)])
set(gca,'YDir','reverse')
xlabel('Mean Absorption Enhancement','Interpreter','Latex')
ylabel('Depth (nm)','Interpreter','Latex')
set(gca,'Fontsize',16)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w');
box on

%with QD layers
d_max = sum(d);
d_temp=0;
for i=1:height(QD_table)
  if d_temp<QD_table(i,2)
      n(2*i-1,:) = retindice_semicond(lambda*1e6,40);
      d(2*i-1)   = (QD_table(i,2)-d_temp-1);
      q(2*i-1)   = 0;
      n(2*i,:)   = retindice_semicond(lambda*1e6,81.1);
      d(2*i)     = 1;
      q(2*i)     = 1;
      d_temp     = QD_table(i,2);
  end
end
n(width(d)+1,:) = retindice_semicond(lambda*1e6,40);
d(width(d)+1)   = (d_max-d_temp);
q(width(d)+1)   = 0;

[~, I_z, ~]=TMM(nb_lambda,lambda,cal_abs,cal_field,transform_E,n_0,n_end,n,d);

I_mean=0;
for j=1:length(q)
    if q(j)>0
        I_mean=I_mean+0.1*I_z{j};
    end
end

figure
plot(lambda*1e6,I_mean,'Linewidth',3);
xlim([min(lambda*1e6) max(lambda*1e6)])
xlabel('$\lambda \: \mathrm{(\mu m)}$','Interpreter','Latex')
ylabel('Mean Absorption Enhancement','Interpreter','Latex')
set(gca,'Fontsize',16)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w');
box on

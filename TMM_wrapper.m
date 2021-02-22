clear all

%% parameters
%wavelength range
lambda=(0.920:0.001:1.340); %Wavelength vector
nb_lambda=width(lambda); %Number of wavelength

%flags
cal_abs=1; % calculate the absorption and field intensity in each layer
cal_field=1; % calculate the field intensity in each layer
transform_E=0; % Transform the data into a function of energy

% structure parameters
back_b = 50;
surf_b = 25;
spacer = 25;
intensity = 12;
height_QD = 1;

%% without QD layers
% List of layers
n_0=ones(1,nb_lambda); %Incident medium
n_end=retindice_metal(lambda,1.9); %Last medium

% list of indices of the stack
n(1,:)=retindice_semicond(lambda,40);

% list of thicknesses of the stack
d(1)=1800;

% calculate normalized field intensity
[I, ~, stack]=TMM(nb_lambda,lambda,cal_abs,cal_field,transform_E,n_0,n_end,n,d);

% I for evaluation
I_mean = I(:,abs(lambda-1.192) < 0.0001);    %I_mean = mean(I,2); % mean I for each hight

% sort in descending order of I 
I_table = cat(2, I_mean, stack');
I_sorted = sortrows(I_table(back_b:end-surf_b,:),'descend');

% set QD layers
QD_table = I_sorted(1,:);
for i=2:height(I_sorted)
    if abs(I_sorted(i,2)-QD_table(:,2))>spacer & I_sorted(i,1)>=intensity
        QD_table = cat(1, QD_table, I_sorted(i,:));
    end
end
QD_table = sortrows(QD_table,2);

%Plot Absorption Enhancement v.s. Depth
drawGraph('Absorption_Enhancement_Depth',I_table,QD_table,d)

%% with QD layers
% List of layers
d_max = sum(d);
d_temp=0;
for i=1:height(QD_table)
  if d_temp<QD_table(i,2)
      n(2*i-1,:) = retindice_semicond(lambda,40);
      d(2*i-1)   = (QD_table(i,2)-d_temp-height_QD);
      q(2*i-1)   = 0;
      n(2*i,:)   = retindice_semicond(lambda,81.1);
      d(2*i)     = height_QD;
      q(2*i)     = 1;
      d_temp     = QD_table(i,2);
  end
end
n(width(d)+1,:) = retindice_semicond(lambda,40);
d(width(d)+1)   = (d_max-d_temp);
q(width(d)+1)   = 0;

% calculate normalized field intensity
[~, I_z, ~]=TMM(nb_lambda,lambda,cal_abs,cal_field,transform_E,n_0,n_end,n,d);

% calculate mean absorption of QD layers
I_mean=0;
for j=1:length(q)
    if q(j)>0
        I_mean=I_mean+I_z{j};
    end
end
I_mean = I_mean/length(find(q==1));

%Plot Absorption Enhancement v.s. Lambda
drawGraph('Absorption_Enhancement_Lambda',lambda,I_mean)


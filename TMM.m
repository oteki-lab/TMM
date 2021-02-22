function [I, I_z, stack]=TMM(nb_lambda,lambda,cal_abs,cal_field,transform_E,n_0,n_end,n,d)
%% Define parameters
h=6.626e-34; %Planck's constant
c=2.998e8; %Speed of light
q=1.602e-19; % Elementary charge
length_step=1e-9; %discretization step of the stack, for intensity calculation

%% Initialization
nb_layers=size(n,1);
length_stack=sum(d);
nb_steps=round(d);
stack=linspace(1,length_stack,sum(nb_steps));
n_tot=[n_0; n; n_end]; %All indices, including media before and after the stack

Omega=cell(1,nb_lambda);
r=zeros(1,nb_lambda);
t=zeros(1,nb_lambda);
thetaz=zeros(nb_layers,sum(nb_steps),nb_lambda);
Delta=cell(nb_layers+1,nb_lambda);
Upsilon=cell(nb_layers,nb_lambda);
I_z=cell(1,nb_layers);
I_poynting=zeros(nb_layers+1,nb_lambda);
Omega_m_prime=cell(nb_layers+1,nb_lambda);
Omega_m=cell(nb_layers+1,nb_lambda);

%% Calculation of R and T
if isempty(n)
    r=(n_0-n_end)./(n_0+n_end);
    t=2*n_0./(n_0+n_end);
    R=abs(r).^2;
    T=(real(n_end)./real(n_0)).*abs(t).^2; % Transmittance
    A=1-R-T; % Absorption of the whole stack
else
    theta=2*pi*n.*d'.*length_step./(lambda/1e6); %Phase shift in a layer
    
    for j=1:nb_lambda
        for i=nb_layers:-1:1
            Delta{i+1,j}=(1/(2*n_tot(i+1,j)))*[n_tot(i+1,j)+n_tot(i+2,j) n_tot(i+1,j)-n_tot(i+2,j); n_tot(i+1,j)-n_tot(i+2,j) n_tot(i+1,j)+n_tot(i+2,j)];
            
            if i==nb_layers
                Omega_m{end,j}=(Delta{end,j}); % transfer matrix between layers 1 and i
            else
                Omega_m{i+1,j}=Delta{i+1,j}*Omega_m_prime{i+1,j};
            end
            Upsilon{i,j}=[exp(-1i*theta(i,j)) 0; 0 exp(1i*theta(i,j))];
            Omega_m_prime{i,j}=Upsilon{i,j}*Omega_m{i+1,j};
        end
        Delta{1,j}=(1/(2*n_tot(1,j)))*[n_tot(1,j)+n_tot(2,j) n_tot(1,j)-n_tot(2,j); n_tot(1,j)-n_tot(2,j) n_tot(1,j)+n_tot(2,j)];
        Omega_m{1,j}=Delta{1,j}*Omega_m_prime{1,j};
        Omega_m_prime{end,j}=[1 0;0 0];
        Omega{j}=Omega_m{1,j}; % Transfer matrix of the system
        r(j)=Omega{j}(2,1)/Omega{j}(1,1); % Reflection coefficient
        t(j)=1/Omega{j}(1,1); % Transmission coefficient
    end
    
    R=abs(r).^2; % Reflectance
    T=(real(n_end)./real(n_0)).*abs(t).^2; % Transmittance
    A=1-R-T; % Absorption of the whole stack
end

%% Calculation of absorptivity for each layer
if cal_abs==1
    for j=1:nb_lambda
        for i=1:nb_layers+1
            I_poynting(i,j)=abs(t(j))^2/real(n_0(j))*real(n_tot(i+1,j)*conj(Omega_m_prime{i,j}(1,1)+Omega_m_prime{i,j}(2,1))*(Omega_m_prime{i,j}(1,1)-Omega_m_prime{i,j}(2,1)));
        end
    end
    Abs_layer=I_poynting(1:nb_layers,:)-I_poynting(2:nb_layers+1,:);
    
    % Getting the energy vector
    if transform_E
        E_A=fliplr(h*c./(q*(lambda/1e6)));
        A_GaAs_E=fliplr(Abs_layer(2,:));
        Abs_layer_E=fliplr(Abs_layer);
        A_total_E=fliplr(A+T);

        figure
        plot(E_A,A_total_E,'linewidth',3)
        xlim([1.3 3.5])
            save('Abs_TMM_AlGaAs_GaAs_HCSC_ELO_MBE2_Q3_2_2_wav_300_1000_habs_20_h_both_AlGaAs_82_Au_mirror.mat','E_A','A_GaAs_E','Abs_layer_E','A_total_E')
    end

    drawGraph('absorption', lambda,A,T) % Plots
end

%% Calculation of intensity of electric field with depth
if cal_field==1
    for j=1:nb_lambda
        for i=1:nb_layers
            for k=1:nb_steps(i)
                thetaz(i,k,j)=(k-1/2)*length_step*2*pi*n(i,j)/(lambda(j)/1e6);
                I_z{i}(k,j)=real(n(i,j))/real(n_0(j))*abs(t(j)*(Omega_m_prime{i,j}(1,1)*exp(1i*thetaz(i,k,j))+Omega_m_prime{i,j}(2,1)*exp(-1i*thetaz(i,k,j))))^2;
                A_z{i}(k,j)=(4*pi*imag(n(i,j))/(lambda(j)/1e6))*I_z{i}(k,j);
                LDOS_z{i}(k,j)=abs(1+Omega_m_prime{i,j}(2,1)*exp(-1i*thetaz(i,k,j))/(Omega_m_prime{i,j}(1,1)*exp(1i*thetaz(i,k,j))))^2;
            end
            for k=1:nb_steps(i)+1
                thetaz_poynting{i}(k,j)=(k-1)*length_step*2*pi*n(i,j)/(lambda(j)/1e6);
                I_poynting_z{i}(k,j)=abs(t(j))^2/real(n_0(j))*real(n_tot(i+1,j)*conj(Omega_m_prime{i,j}(1,1)*exp(1i*thetaz_poynting{i}(k,j))+Omega_m_prime{i,j}(2,1)*exp(-1i*thetaz_poynting{i}(k,j)))*(Omega_m_prime{i,j}(1,1)*exp(1i*thetaz_poynting{i}(k,j))-Omega_m_prime{i,j}(2,1)*exp(-1i*thetaz_poynting{i}(k,j))));
            end
            for k=1:nb_steps(i)
                A_poynting_z{i}(k,j)=(I_poynting_z{i}(k,j)-I_poynting_z{i}(k+1,j))/length_step;
            end
        end
    end

    I=I_z{1};
    A_local=A_z{1};
    LDOS=LDOS_z{1};
    A_poynting_local=A_poynting_z{1};
    I_poynting_local=I_poynting_z{1};
    for i=2:nb_layers
        I=cat(1,I,I_z{i});
        A_local=cat(1,A_local,A_z{i});
        LDOS=cat(1,LDOS,LDOS_z{i});
        A_poynting_local=cat(1,A_poynting_local,A_poynting_z{i});
        I_poynting_local=cat(1,I_poynting_local,I_poynting_z{i}(2:end,:));
    end
    
    I_poynting_local=(I_poynting_local(2:end,:)+I_poynting_local(1:end-1,:))/2;
    
    drawGraph('normalized_field_intensity', lambda,stack,I,nb_layers,d) % Plots
end
function varargout = drawGraph(request, varargin)
  switch request
     case 'normalized_field_intensity'
       [varargout{1:nargout}] = normalized_field_intensity(varargin{:});
     case 'absorption'
       [varargout{1:nargout}] = absorption(varargin{:});
     case 'Absorption_Enhancement_Depth'
       [varargout{1:nargout}] = Absorption_Enhancement_Depth(varargin{:});
     case 'Absorption_Enhancement_Lambda'
       [varargout{1:nargout}] = Absorption_Enhancement_Lambda(varargin{:});
  end
end
  
 function normalized_field_intensity(lambda,stack,I,nb_layers,d)
    figure
    contourf(lambda,stack/1e3,I);
    hold all
    for i=1:nb_layers-1
        plot([lambda(1),lambda(end)],[sum(d(1:i)) sum(d(1:i))]/1e3,'linewidth',1,'color','w')
    end
    xlabel('$\lambda \: \mathrm{(\mu m)}$','Interpreter','Latex')
    ylabel('depth $\: \mathrm{(\mu m)}$','Interpreter','Latex')
    set(gca,'Fontsize',18)
    set(gca,'YDir','reverse')
    set(gca,'XMinorTick','on','YMinorTick','on')
%     set(gca,'XAxisLocation','top')
%     set(gca,'xticklabel',[])
    set(gcf,'color','w');
    colormap('parula')
%     h.LevelList=[0 1e-3 1e-2 1e-1 0.25 0.5 1 2 4 8];
    box on
    c=colorbar;
    c.Label.String = 'Normalized field intensity';
    ylim([0 inf])
 end

 function absorption(lambda,A,T)
    figure
    plot(lambda,A+T,'linewidth',3);
    % plot(lambda,1-R);
    hold on
    plot(lambda,A,'linewidth',2);
    plot(lambda,T,'linewidth',2);
    ylim([0 1])
    xlim([min(lambda) max(lambda)])
    xlabel('$\lambda \: \mathrm{(\mu m)}$','Interpreter','Latex')
    ylabel('$A$','Interpreter','Latex')
    set(gca,'Fontsize',16)
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gcf,'color','w');
    legend({'Total','QDs','Ag'})
    box on
 end

 function Absorption_Enhancement_Depth(I_table,QD_table,d)
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
 end

function Absorption_Enhancement_Lambda(lambda,I_mean)
    figure
    plot(lambda,I_mean,'Linewidth',3);
    xlim([min(lambda) max(lambda)])
    xlabel('$\lambda \: \mathrm{(\mu m)}$','Interpreter','Latex')
    ylabel('Mean Absorption Enhancement','Interpreter','Latex')
    set(gca,'Fontsize',16)
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gcf,'color','w');
    box on
end
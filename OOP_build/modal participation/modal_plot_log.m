function mode_fig=modal_plot_log(phi_jkn,all_harmonics,state_number,state_name,mode_names,modes_to_plot,colors_to_plot)

fontsize_label=11;
mode_fig = figure;
t=tiledlayout(length(modes_to_plot),1);
xlabel(t,'Harmonic','Interpreter','latex')
ylabel(t,'$-log_{10}(1-\phi)$','Interpreter','latex')
title(t,state_name,'Interpreter','latex')
t.TileSpacing='compact';

% plot of all modes indicated in modes_to_plot

nexttile
    for i=1:size(phi_jkn,3)
        phi_to_plot(i) = -log10(1-phi_jkn(state_number,modes_to_plot(1),i));
    end
    stem(all_harmonics,phi_to_plot,'o','MarkerSize',6,'LineWidth',2,'Color',colors_to_plot(1,:))
    xticks((all_harmonics(1)-1:all_harmonics(end)+1))
    %yticks((-2:0.5:0))
    grid on
    xlim([all_harmonics(1)-1  all_harmonics(end)+1])
     %ylim([0 1])
    ylabel(mode_names(1),'Interpreter','latex','FontSize',fontsize_label)
    legend('LTP','interpreter','latex','Location','northwest')

for m=2:length(modes_to_plot)
    % plot of a single mode
    nexttile
    for i=1:size(phi_jkn,3)
        phi_to_plot(i) = -log10(1-phi_jkn(state_number,modes_to_plot(m),i));
    end
    stem(all_harmonics,phi_to_plot,'o','MarkerSize',6,'LineWidth',2,'Color',colors_to_plot(m,:))
    xticks((all_harmonics(1)-1:all_harmonics(end)+1))
    %yticks((0:0.25:1))
    grid on
    xlim([all_harmonics(1)-1  all_harmonics(end)+1])
    %ylim([0 1])
    ylabel(mode_names(m),'Interpreter','latex','FontSize',fontsize_label)
    
end
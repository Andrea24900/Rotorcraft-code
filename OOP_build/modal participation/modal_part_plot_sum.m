function [modal_plot,phi_1mode1state] = modal_part_plot_sum (A_handle_all,char_exp,eigenvectors,...
    rotor_OMEGA_vector,Nh,mode,state,figure_number,ylabel_text,nolegend)
colors

for i =1: length(rotor_OMEGA_vector)

    phi_1mode1state (:,i) = modal_participation_par(...
    A_handle_all,char_exp(:,i),eigenvectors(:,:,i),Nh,mode,state,rotor_OMEGA_vector(i));
    k=1;
    for j=1:ceil(height(phi_1mode1state)/2)
        phi_1mode1state_sum (k,i) = phi_1mode1state(j,i) + phi_1mode1state(1+height(phi_1mode1state)-j,i);
        k= k+1;
    end
end


color_vect = [blue_mat; orange_mat; red_mat; green_mat; [0 0 0]; [0 1 0]; purple_mat; [1 0 0];[0 0 1]];
modal_plot = figure(figure_number);

set(gcf,'Units','normalized','Position',[0.2, 0.2, 0.3, 0.4])
patch([210 305 305 210],[0 0 1 1],'r','facealpha',0.05)
hold on
grid on
plot (rotor_OMEGA_vector*60/2/pi,phi_1mode1state_sum(1,:),'LineWidth',line_width,'Color',color_vect(1,:))

for i = 2:height(phi_1mode1state_sum)-1
    plot (rotor_OMEGA_vector*60/2/pi,phi_1mode1state_sum(i,:),'-','LineWidth',line_width,'Color',color_vect(i,:))
end
i = height(phi_1mode1state_sum);
    plot (rotor_OMEGA_vector*60/2/pi,phi_1mode1state_sum(i,:),'-.','LineWidth',line_width,'Color',color_vect(i,:))

if nolegend~=1
legend({'Instability','$n=\pm 4$','$n=\pm 3$','$n=\pm 2$','$n=\pm 1$','$n=0$'},'Interpreter',...
    'latex','FontSize',fontsize_legend,'location','eastoutside');
end
xlim([0 400])
ylim([0 1])
axis([0 400 0 1])
set(gca,'Position',[0.1, 0.15, 0.5, 0.7])
%title(title_text,'Interpreter','latex','FontSize',14)
xlabel('$\Omega$ [rpm]','Interpreter','latex','FontSize',fontsize_label)
ylabel(ylabel_text,'Interpreter','latex','FontSize',18)
end
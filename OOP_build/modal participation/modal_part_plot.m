function [modal_plot,phi_1mode1state] = modal_part_plot (A_handle_all,char_exp,eigenvectors,...
    rotor_OMEGA_vector,Nh,mode,state,figure_number,ylabel_text,nolegend)
colors

for i =1: length(rotor_OMEGA_vector)
phi_1mode1state (:,i) = modal_participation_par(...
    A_handle_all,char_exp(:,i),eigenvectors(:,:,i),Nh,mode,state,rotor_OMEGA_vector(i));
end

color_vect = [blue_mat; orange_mat; red_mat; green_mat; [0 0 0]; [0 1 0]; purple_mat; [1 0 0];[0 0 1]];
modal_plot = figure(figure_number);

set(gcf,'Units','normalized','Position',[0.2, 0.2, 0.3, 0.4])
patch([210 305 305 210],[0 0 1 1],'r','facealpha',0.05)
hold on
grid on
plot (rotor_OMEGA_vector*60/2/pi,phi_1mode1state(1,:),'LineWidth',line_width,'Color',color_vect(1,:))

for i = 2:floor(height(phi_1mode1state)/2)
    plot (rotor_OMEGA_vector*60/2/pi,phi_1mode1state(i,:),'-','LineWidth',line_width,'Color',color_vect(i,:))
end
i = ceil(height(phi_1mode1state)/2);
    plot (rotor_OMEGA_vector*60/2/pi,phi_1mode1state(i,:),'-.','LineWidth',line_width,'Color',color_vect(i,:))

for i = ceil(height(phi_1mode1state)/2)+1:height(phi_1mode1state)
    plot (rotor_OMEGA_vector*60/2/pi,phi_1mode1state(i,:),'--','LineWidth',line_width,'Color',color_vect(i,:))
end
if nolegend~=1
legend({'Instability','$n=-4$','$n=-3$','$n=-2$','$n=-1$','$n=0$','$n=+1$','$n=+2$','$n=+3$','$n=+4$'},'Interpreter',...
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
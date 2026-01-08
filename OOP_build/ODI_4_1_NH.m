clear all
close all
clc
%set 3 harmonics and HD to see the effect of many harmonics in the
%decomposition


rotor_build = rotor_build.build_all();

modes = stability_analysis.run_stability(rotor_build);
%%

%modal_participation = modal_participation_analysis(modes);

%%
close all

%my_plot.plot_mod_part(modal_participation.modal_participation,8,4,...
 %   {'$\phi_{Z_{1c},4,n}$'},'XLimits',[50 400],'PlotSum',false)


modes_to_plot_0 = [9 10 11 12 5 6 3 4 7 8 1 2];

handle=my_plot.plot_damping_generic(modes.modal_solution,'Xlimits',[0 400]);
%my_plot.plot_damping_order(modes.modal_solution,'Modes',modes_to_plot_0,'FigureHandle',handle)
my_plot.plot_damping_order(modes.modal_solution,'Modes',modes_to_plot_0)
% %% first harmonic (modes from 13 to 36)
% modes_to_plot_1c = [15 16 13 14 19 20 23 24 17 18 21 22];
%
% my_plot.plot_damping_order(modes.modal_solution,'Modes',modes_to_plot_1c)
% my_plot.plot_mod_part(modal_participation.modal_participation,8,9,...
%     {'$\phi_{Z_{1c},9,n}$'},'XLimits',[50 400],'PlotSum',false)
% %handle_2=my_plot.plot_frequency_generic(modes.modal_solution,'Xlimits',[35 200],'Ylimits',[0 45]);
% %my_plot.plot_frequency_order(modes.modal_solution,'Modes',(1:24),'FigureHandle',handle_2)
%
% modes_to_plot_1s = 24*ones(1,12)+modes_to_plot_0;
%
%
% my_plot.plot_damping_order(modes.modal_solution,'Modes',modes_to_plot_1s)
% my_plot.plot_mod_part(modal_participation.modal_participation,8,16,...
%     {'$\phi_{Z_{1c},16,n}$'},'XLimits',[50 400],'PlotSum',false)

clear all 
close all
clc



rotor_build = rotor_build.build_all();

modes = stability_analysis.run_stability(rotor_build);
%%

modal_participation = modal_participation_analysis(modes);

%%
close all

my_plot.plot_mod_part(modal_participation.modal_participation,8,4,...
    {'$\phi_{Z_{1c},4,n}$'},'XLimits',[50 400],'PlotSum',false)


modes_to_plot = (1:12);

%handle=my_plot.plot_damping_generic(modes.modal_solution,'Xlimits',[0 400]);
%my_plot.plot_damping_order(modes.modal_solution,'Modes',modes_to_plot,'FigureHandle',handle)
my_plot.plot_damping_order(modes.modal_solution,'Modes',modes_to_plot)
  
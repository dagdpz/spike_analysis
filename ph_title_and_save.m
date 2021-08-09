function ph_title_and_save(figure_handle,filename,plot_title,keys)
mtit(figure_handle,  plot_title, 'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none');
%stampit;

wanted_size=[50 30];
set(figure_handle, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
if keys.plot.export
    export_fig([keys.path_to_save, filesep, filename], '-pdf','-transparent') % pdf by run
    close(figure_handle);
end
end
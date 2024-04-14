function save_all_figs(folder_name) 

    % folder_name = 'outputs';   % Your destination folder
    fig_list    = findobj(allchild(0), 'flat', 'Type', 'figure');

    % loop through all figs and save 
    for i_fig = 1:length(fig_list)
      fig_handle = fig_list(i_fig);
      ax_handle  = fig_handle.CurrentAxes; 
      fig_name   = get( get( ax_handle, 'title' ) , 'string');
      fig_save   = fullfile(folder_name, [ fig_name, '.png' ] ) ; 
      saveas(fig_handle, fig_save );
    end

end 
function Record_RotatingGraph(figure_handle,save_name)
    arguments
        figure_handle matlab.ui.Figure
        save_name
    end
    x = pi/4: pi/180*0.1 : pi/2 + pi/4;

    frames(1000) = struct('cdata', [], 'colormap', []);
    frames(1) = getframe(figure_handle);
    v = VideoWriter(save_name,'Motion JPEG AVI');
    v.Quality = 100;
    open(v);
    for k = 1 : length(x)
        view([20+0.1*k,20])
        % pause(0);
        frames(k+1) = getframe(figure_handle);
    end
    writeVideo(v,frames(1:(k-1)));
    close(v);
end
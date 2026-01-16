function plot_movie(name, usol, curve, ps)
    % Function to plot and save a movie from the u solution
    % Inputs:
    %   name : File name for the movie. Should end in .avi
    %   usol : Solution to the time domain problem
    %   curve: The curve scattered off of
    %   ps   : The parameters.
    lim = max(abs(real(usol(:,:,:))), [],"all");
    % lim = 0.5;

    v = VideoWriter(name);
    open(v);
    figure(3)
    clf

    frameArray = repmat(getframe(figure(3)), ps.numt, 1);
    figure(3)
    

   
    offset = 0;
    for t=1:ps.numt - offset
        clf
        colorbar
        ax = gca;
        % ax.CLim = [-lim, lim];
        ax.CLim = [0,lim];
        hold on
        % hold on
        % We transpose here because x in in the rows, and y in the columns
        toplot = abs(real(usol(:,:,t+offset))).';
        imagesc(ps.xlims,ps.ylims,toplot,'Interpolation','bilinear')
        % pcolor(ps.xs,ps.ys,real(toplot))
        if ps.is_open_curve
            plot(curve.X,curve.Y,'k', 'Linewidth',2);
        else
            curve.draw()
        end
        shading interp;
        colormap(jet)
        title("t = " + num2str(ps.ts(t+offset)));
        xlim(ps.xlims);
        ylim(ps.ylims);

        
        
        if ps.is_open_curve
            plot(curve.X, curve.Y, 'k')
        else
            plot(curve.xs(1,:),curve.xs(2,:),'k')
        end
        
        frameArray(t) = getframe(gcf);
        % hold off
    end
    hold off
    

    tic
    writeVideo(v,frameArray)
    wtime = toc;

    close(v);

end

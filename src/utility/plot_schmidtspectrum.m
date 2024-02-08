function [ax] = plot_schmidtspectrum(S, ax, kwargs)

arguments
    S
    ax = []
    kwargs.SymmetrySort = true
    kwargs.ExpandQdim = false
    kwargs.Marker = '.'
end

    if isempty(ax)
        figure;
        ax = gca;
    end

    if kwargs.Marker == '.'
        kwargs.MarkerSize = 10;
    else
        kwargs.MarkerSize = 6;
    end
    
    %hold(ax, 'off');
    [svals, charges] = matrixblocks(S);
    if kwargs.ExpandQdim
        for i = 1:length(svals)
            svals{i} = reshape(repmat(svals{i}, 1, qdim(charges(i))), [], 1);
        end
    end
    ctr = 0;
    labels = arrayfun(@string, charges, 'UniformOutput', false);
    lengths = cellfun(@length, svals);
    ticks = cumsum(lengths);
    if kwargs.SymmetrySort
        for i = 1:length(svals)
            semilogy(ax, ctr+(1:lengths(i)), svals{i}, 'Marker', kwargs.Marker, 'MarkerSize', kwargs.MarkerSize, 'Color', colors(i));
            if i == 1, hold(ax, 'on'); end
            ctr = ctr + lengths(i);
        end
        set(ax, 'Xtick', ticks, 'fontsize', 10, ...
            'XtickLabelRotation', 60, 'Xgrid', 'on');
    else
        [~, p] = sort(vertcat(svals{:}), 'descend');
        p = invperm(p);
        for i = 1:length(svals)
            semilogy(ax, p(ctr+(1:lengths(i))), svals{i}, 'Marker', kwargs.Marker, 'MarkerSize', kwargs.MarkerSize, 'Color', colors(i));
            if i == 1, hold(ax, 'on'); end
            ctr = ctr + lengths(i);
        end
        
    end
    legend(ax, arrayfun(@(x) line([0,1],[0,0],'Color', colors(x), 'Marker', '.', 'MarkerSize', 10, 'linestyle', 'none'),1:numel(labels)), labels, 'Location', 'Best')
    set(ax, 'TickLabelInterpreter', 'latex');
    xlim(ax, [1 - 1e-8 ticks(end) + 1e-8]);
    hold off
    
    linkaxes(ax, 'y');

end


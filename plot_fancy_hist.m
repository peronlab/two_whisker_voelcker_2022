function plot_fancy_hist (ax, dat_vec, bins, colr, lw)
    [n x] = hist(dat_vec, bins);
    n=n/nansum(n);

    plot_x = reshape(repmat(x,2,1),[],1);
    plot_y = plot_x*nan;
    n = [nan n];
    for i=2:2:length(plot_y)
        ii = floor(i/2);
        plot_y([i-1 i]) = [n(ii) n(ii+1)];
    end
    plot(ax, plot_x, plot_y, '-', 'Color', colr, 'LineWidth', lw);


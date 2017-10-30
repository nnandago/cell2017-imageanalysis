function background = find_background(centroid, im, segim, bin)
    global imx imy
    [upperbound_x, upperbound_y, ~, ~] = find_bounds([centroid(2) - 30/bin, centroid(1) - 30/bin], [60/bin, 60/bin], [imx imy]);
    neighborhood_im = im(uint16(max(1, centroid(2) - 30/bin)):upperbound_x, uint16(max(1, centroid(1) - 30/bin)):upperbound_y);
    neighborhood_segim = segim(uint16(max(1, centroid(2) - 30/bin)):upperbound_x, uint16(max(1, centroid(1) - 30/bin)):upperbound_y);
    background = median(double(neighborhood_im(neighborhood_segim == 0)));
end
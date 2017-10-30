function sub_im0 = find_subimage(im, bound, mask)
    start_row = uint16(bound(2)); start_col = uint16(bound(1));
    sub_im0 = im(start_row: start_row + size(mask, 1) - 1, start_col: start_col + size(mask, 2) - 1).*uint16(mask);
end

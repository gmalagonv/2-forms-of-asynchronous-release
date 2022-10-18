function [Norm_Fac, ectopic] = func_Norm_Fac(Asyn_loc_XY, centerPhyXY, area_point)
[crossx, crossy] = polyxpoly([area_point(:, 1); area_point(1, 1)], [area_point(:, 2); area_point(1, 2)], ...
    [centerPhyXY(1, 1), Asyn_loc_XY(1, 1)], [centerPhyXY(1, 2), Asyn_loc_XY(1, 2)]);
expan_edge = 25;                                                                                                                % <======= addtional boundary 25nm
if isempty([crossx, crossy]) || pdist([Asyn_loc_XY; centerPhyXY]) - pdist([crossx, crossy; centerPhyXY]) < expan_edge
    temp_Asyn_loc_XY = Asyn_loc_XY;
    while isempty([crossx, crossy]) || pdist([temp_Asyn_loc_XY; centerPhyXY]) - pdist([crossx, crossy; centerPhyXY]) < expan_edge
        temp_Asyn_loc_XY = (temp_Asyn_loc_XY - centerPhyXY) * 2 + centerPhyXY;
        [crossx, crossy] = polyxpoly([area_point(:, 1); area_point(1, 1)], [area_point(:, 2); area_point(1, 2)], ...
            [centerPhyXY(1, 1), temp_Asyn_loc_XY(1,1)], [centerPhyXY(1,2), temp_Asyn_loc_XY(1, 2)]);                                 % new temporal cross x and y
    end
    Norm_Fac = pdist([crossx, crossy; centerPhyXY]) + expan_edge;
    ectopic = false;
else
    ectopic = true;
    [crossx, crossy] = polyxpoly([area_point(:, 1); area_point(1, 1)], [area_point(:, 2); area_point(1, 2)], ...
        [centerPhyXY(1, 1), Asyn_loc_XY(1,1)], [centerPhyXY(1,2), Asyn_loc_XY(1, 2)]);                                 % new temporal cross x and y
    Norm_Fac = pdist([crossx, crossy; centerPhyXY]) + expan_edge;
end
end
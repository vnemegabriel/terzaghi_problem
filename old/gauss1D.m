function [w,gp] = gauss1D(n)

switch n
    case 1
        w = 2;
        gp = 0;
    case 2
        w = [1;1];
        a = sqrt(3)/3;
        gp = [-a;a];
    case 3
        w = [5;8;5]/9;
        a = sqrt(3/5);
        gp = [-a;0;a];

end

end
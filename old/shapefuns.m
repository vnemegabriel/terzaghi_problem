function Ni = shapefuns(locationGaussPoints,eleType)

xi = locationGaussPoints(:,1);
eta = locationGaussPoints(:,2);

switch eleType

    case 'Q4'    
        N1 = 0.25*(1-xi)*(1-eta);
        N2 = 0.25*(1+xi)*(1-eta);     
        N3 = 0.25*(1+xi)*(1+eta);
        N4 = 0.25*(1-xi)*(1+eta);
        Ni = [N1 N2 N3 N4];


    case 'Q8'
        N8 = 1/2 * (1 - xi) * (1 - eta^2);
        N7 = 1/2 * (1 - xi^2) * (1 + eta);
        N6 = 1/2 * (1 + xi) * (1 - eta^2);
        N5 = 1/2 * (1 - xi^2) * (1 - eta);
        N4 = 1/4  * (1 - xi) * (1 + eta) - 1/2 * (N7 + N8);
        N3 = 1/4  * (1 + xi) * (1 + eta) - 1/2 * (N6 + N7);
        N2 = 1/4  * (1 + xi) * (1 - eta) - 1/2 * (N5 + N6);
        N1 = 1/4  * (1 - xi) * (1 - eta) - 1/2 * (N8 + N5);

        % Ni = [N1 N5 N2 N6 N3 N7 N4 N8];

        Ni = [N1 N2 N3 N4 N5 N6 N7 N8];

    case 'CST'
        
        N1 = xi;
        N2 = eta;
        N3 = 1 - xi - eta;
        Ni = [N1 N2 N3];

    
end


end


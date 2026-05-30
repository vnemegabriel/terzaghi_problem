function dN = shapefunsDer(locationGaussPoints,eleType)

xi = locationGaussPoints(1);
eta = locationGaussPoints(2);
% Solo son un elemento cada uno, no son arrays

switch eleType

    case 'Q4'        
        %xi derivatives 
        dN1_dxi = -0.25*(1-eta);
        dN2_dxi = 0.25*(1-eta);  
        dN3_dxi = 0.25*(1+eta);  
        dN4_dxi = -0.25*(1+eta); 
        
        %eta derivatives
        dN1_deta = -0.25*(1-xi);
        dN2_deta = -0.25*(1+xi);       
        dN3_deta = 0.25*(1+xi);
        dN4_deta = 0.25*(1-xi);

        dN_xi = [dN1_dxi dN2_dxi dN3_dxi dN4_dxi];
        dN_eta = [dN1_deta dN2_deta dN3_deta dN4_deta];
        
        dN = [dN_xi;dN_eta];
          

    case 'Q8'
        %xi derivatives
        dN8_dxi = 0.5 * (-1) * (1 - eta^2);
        dN7_dxi = 0.5 * (-2*xi) * (1 + eta);
        dN6_dxi = 0.5 * (1) * (1 - eta^2);
        dN5_dxi = 0.5 * (-2*xi) * (1 - eta);
        dN4_dxi = 0.25 * (-1) * (1 + eta) - 0.5 * (dN7_dxi + dN8_dxi);
        dN3_dxi = 0.25 * (1) * (1 + eta) - 0.5 * (dN6_dxi + dN7_dxi);
        dN2_dxi = 0.25 * (1) * (1 - eta) - 0.5 * (dN5_dxi + dN6_dxi);
        dN1_dxi = 0.25 * (-1) * (1 - eta) - 0.5 * (dN5_dxi + dN8_dxi);

        %eta derivatives
        dN8_deta = 0.5 * (-2*eta) * (1 - xi);
        dN7_deta = 0.5 * (1) * (1 - xi^2);
        dN6_deta = 0.5 * (-2*eta) * (1 + xi);
        dN5_deta = 0.5 * (-1) * (1 - xi^2);
        dN4_deta = 0.25 * (1) * (1 - xi) - 0.5 * (dN7_deta + dN8_deta);
        dN3_deta = 0.25 * (1) * (1 + xi) - 0.5 * (dN6_deta + dN7_deta);
        dN2_deta = 0.25 * (-1) * (1 + xi) - 0.5 * (dN5_deta + dN6_deta);
        dN1_deta = 0.25 * (-1) * (1 - xi) - 0.5 * (dN5_deta + dN8_deta);

        dN_xi = [dN1_dxi dN2_dxi dN3_dxi dN4_dxi dN5_dxi dN6_dxi dN7_dxi dN8_dxi];
        dN_eta = [dN1_deta dN2_deta dN3_deta dN4_deta dN5_deta dN6_deta dN7_deta dN8_deta];

        dN = [dN_xi;dN_eta];

    case 'CST'
        
        % xi derivatives
        dN1_dxi = 1;
        dN2_dxi = 0;
        dN3_dxi = -1;
        
        % eta derivatives
        dN1_deta = 0;
        dN2_deta = 1;
        dN3_deta = -1;
        
        dN_xi = [dN1_dxi dN2_dxi dN3_dxi];
        dN_eta = [dN1_deta dN2_deta dN3_deta];

        dN = [dN_xi; dN_eta]; % (2x3 matrix)
end

end
function K = getStiffnessMatrix(nodesElem,eleType,t,E,nPointsGauss) 

npg_xi = nPointsGauss(1);
npg_eta = nPointsGauss(2);

[weights_xi,lgp_xi] = gauss1D(npg_xi);
[weights_eta,lgp_eta] = gauss1D(npg_eta);

nDofNod = 2; %2D
nNodEle = size(nodesElem,1);
nDofEle = nDofNod*nNodEle;

switch eleType

    case 'Q4'
        K = zeros(nDofEle);

        for ipg_xi = 1:npg_xi
            for ipg_eta = 1:npg_eta

            xi = lgp_xi(ipg_xi);
            eta = lgp_eta(ipg_eta);

            dN_iso = shapefunsDer([xi eta],'Q4');
            
            Jacobian = dN_iso*nodesElem;

            assert(det(Jacobian)>0)

            dNxy = Jacobian\dN_iso;
            
            B = zeros(3,nDofEle);
            B(1, 1:2:end) = dNxy(1,:);
            B(2, 2:2:end) = dNxy(2,:);
            B(3, 1:2:end) = dNxy(2,:);
            B(3, 2:2:end) = dNxy(1,:);

            %Integrand
            K = K + weights_xi(ipg_xi)*weights_eta(ipg_eta)*t*B'*E*B*det(Jacobian); 
            end
        end
    case 'Q8'
        K = zeros(nDofEle);

        for ipg_xi = 1:npg_xi
            for ipg_eta = 1:npg_eta

            xi = lgp_xi(ipg_xi);
            eta = lgp_eta(ipg_eta);

            dN_iso = shapefunsDer([xi eta],'Q8');
            
            Jacobian = dN_iso*nodesElem;

            assert(det(Jacobian)>0)

            dNxy = Jacobian\dN_iso;

            %Strain displacement B matrix 
            B = zeros(3,nDofEle);
            B(1, 1:2:end) = dNxy(1,:);
            B(2, 2:2:end) = dNxy(2,:);
            B(3, 1:2:end) = dNxy(2,:);
            B(3, 2:2:end) = dNxy(1,:);

            %Integrand
            K = K + weights_xi(ipg_xi)*weights_eta(ipg_eta)*t*B'*E*B*det(Jacobian); 
            end
        end
            
    case 'CST'
        dN_iso = shapefunsDer([0 0],'CST'); 
        
        Jacobian = dN_iso*nodesElem; 
        detJ = det(Jacobian);
        
        Area = detJ / 2;
        
        % Derivatives in global (x,y) coordinates (2x3)
        dNxy = Jacobian\dN_iso; 
        
        B = zeros(3,nDofEle); % nDofEle is 6 for CST
        B(1, 1:2:end) = dNxy(1,:); % dN/dx
        B(2, 2:2:end) = dNxy(2,:); % dN/dy
        B(3, 1:2:end) = dNxy(2,:); % dN/dy
        B(3, 2:2:end) = dNxy(1,:); % dN/dx

        K = B' * E * B * t * Area;
end

end



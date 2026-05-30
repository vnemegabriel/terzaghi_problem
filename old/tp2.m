%% Trabajo practico #2 
clear 
clc
close all
%% Carga de Malla y Parámetros
mesh = lectorMalla('malla_presa_invertida.dat','xy');

nDofNod = 2;
nNod = size(mesh.nodes,1);
nElem = size(mesh.elements,1);
fprintf('Malla cargada: %d Nodos, %d Elementos.\n', nNod, nElem);
%% Propiedades del material
density = 2400;     % [kg/m^3] 
E = 21000;  % [MPa] = [N/mm^2]
v = 0.2;    % Poisson
t = 1;   % [mm] 
C = E / ((1+v)*(1-2*v)) * [1-v  v   0; ...
                           v  1-v   0; ...
                           0   0  (1-2*v)/2]; % Constitutivo (Plane Strain)
%% Parametros
g = 9.81 ;           % [m/s^2]
body_force_y = -density * g * (1/1000)^3; % [N/mm^3]
b_vec = [0; body_force_y]; 
% --- Presión de Agua ---
rho_water = 1000;   % [kg/m^3] - Densidad Agua
pressure_coeff = rho_water * g * (1/1000)^3; % <-- CORRECTO
water_level_meters = 28;
max_y_water = water_level_meters * 1000; % Nivel del agua en [mm]

%% --- 4. Puntos de Gauss ---
npg_vol = 2; % Para integral de área 2D (Fuerza Másica)
[w_gp_vol, lgp_vol] = gauss1D(npg_vol); 
npg_line = 2; % Para integral de línea 1D (Presión)
[w_gp_line, lgp_line] = gauss1D(npg_line);
%% Assembly de K y F
F = zeros(nNod, nDofNod); % Vector de fuerzas (mapeado a filas)
nSparse = nElem * (8*nDofNod)^2; % Estimación generosa para Q8
iRowSparse = zeros(nSparse, 1);
iColSparse = zeros(nSparse, 1);
valueSparseK = zeros(nSparse, 1);
iloc = 0; % Índice para ensamblaje sparse
for k = 1:nElem
    localNodes = mesh.nodes(mesh.elements{k}, :); % Coordenadas (de las filas correctas)
    nNodElem = numel(mesh.elements{k});
    
    elemType = '';
    edges = [];
    if nNodElem == 4
        elemType = 'Q4';
        edges = [1 2; 2 3; 3 4; 4 1]; % Nodos locales [n1, n2] (Corregido)
    elseif nNodElem == 3
        elemType = 'CST';
        edges = [1 2; 2 3; 3 1];
    elseif nNodElem == 8
        elemType = 'Q8';
        edges = [1 2 5; 2 3 6; 3 4 7; 4 1 8]; % Nodos locales [n1, n2, n_mid]
    else
        continue; % Omitir elemento desconocido
    end
    % --- 5.A. Cálculo de Matriz de Rigidez (kEle) ---
    kEle = getStiffnessMatrix(localNodes, elemType, t, C, [npg_vol npg_vol]);
    
    % --- 5.B. Cálculo de Vector de Fuerzas (Fe) ---
    Fe_body = zeros(nNodElem*nDofNod,1);
    Fe_pressure = zeros(nNodElem*nDofNod,1);
    
    % --- Fuerza Másica (Integral de Área 2D) ---
    if strcmp(elemType, 'CST')
        % Para CST, la integral se puede hacer analíticamente (1 punto)
        N_vec = [1/3, 1/3, 1/3]; % Centroide
        dN_iso = shapefunsDer(N_vec(1:2), elemType); % dN en centroide
        Jacobian = dN_iso * localNodes;
        detJ = det(Jacobian)/2; % detJ es 2*Area
        
        N_matrix_body = zeros(2, nNodElem * nDofNod);
        N_matrix_body(1, 1:2:end) = N_vec; % Componentes X
        N_matrix_body(2, 2:2:end) = N_vec; % Componentes Y
        
        Fe_body = (N_matrix_body' * b_vec) * t * detJ; % * 1 (peso de 1 punto)
    else
        % Para Q4 y Q8, usar bucle de Gauss
        for i = 1:npg_vol
            for j = 1:npg_vol
                xi = lgp_vol(i); eta = lgp_vol(j);
                w_i = w_gp_vol(i); w_j = w_gp_vol(j);
                
                N_vec = shapefuns([xi eta], elemType); 
                dN_iso = shapefunsDer([xi eta], elemType); 
                Jacobian = dN_iso * localNodes; 
                detJ = det(Jacobian);
                
                N_matrix_body = zeros(2, nNodElem * nDofNod);
                N_matrix_body(1, 1:2:end) = N_vec; % Componentes X
                N_matrix_body(2, 2:2:end) = N_vec; % Componentes Y
                
                Fe_body = Fe_body + (N_matrix_body' * b_vec) * t * detJ * w_i * w_j;
            end
        end
    end
    
    % Presion del agua
    for edge_idx = 1:size(edges, 1)
        edge_local_nodes_idx = edges(edge_idx, :); % Índices locales (ej: [1 2 5])
        edge_global_nodes_rows = mesh.elements{k}(edge_local_nodes_idx);
        edge_coords = mesh.nodes(edge_global_nodes_rows, :);
        
        % Chequear si el borde está en la cara izquierda (x=0)
        is_on_left_face = all(abs(edge_coords(:, 1)) < 1e-6); % Usar tolerancia
        
        if is_on_left_face
            nNodEdge = numel(edge_local_nodes_idx);
            for gp_idx = 1:npg_line
                zeta = lgp_line(gp_idx); % Coord. 1D de -1 a 1
                w_gp = w_gp_line(gp_idx);
                
                N_edge = []; dN_edge_dzeta = [];
                if nNodEdge == 2 % Borde lineal (Q4, CST)
                    N_edge = [0.5*(1-zeta), 0.5*(1+zeta)];
                    dN_edge_dzeta = [-0.5, 0.5];
                elseif nNodEdge == 3 % Borde cuadrático (Q8)
                    % N_1, N_2, N_mid (en orden [n1, n2, n_mid])
                    N_edge = [-0.5*zeta*(1-zeta), 0.5*zeta*(1+zeta), (1-zeta^2)];
                    dN_edge_dzeta = [-0.5 + zeta, 0.5 + zeta, -2*zeta];
                end
                
                gp_global_coords = N_edge * edge_coords;
                dCoords_dzeta = dN_edge_dzeta * edge_coords;
                J_1D = norm(dCoords_dzeta); % Jacobiano 1D (dl/dzeta)
                
                current_y = gp_global_coords(2);
                depth = max_y_water - current_y;
                
                if depth > 0
                    pressure = pressure_coeff * depth; % [N/mm^2]
                    normal_vec = [1; 0]; % Normal apunta en +X
                    
                    N_matrix_pressure = zeros(2, nNodElem * nDofNod);
                    for n = 1:nNodEdge
                        node_local_idx = edge_local_nodes_idx(n); % Índice (1-8)
                        dof_x = (node_local_idx - 1) * nDofNod + 1;
                        dof_y = (node_local_idx - 1) * nDofNod + 2;
                        
                        N_matrix_pressure(1, dof_x) = N_edge(n);
                        N_matrix_pressure(2, dof_y) = N_edge(n);
                    end
                    
                    Fe_pressure = Fe_pressure + (N_matrix_pressure' * normal_vec) * pressure * J_1D * w_gp * t;
                end
            end
        end
    end    
    % Ensamblaje de K (Sparse)
    for i1 = 1:nNodElem
        for i2 = 1:nDofNod
            iCol = i2+(i1-1)*nDofNod;
            iColGlobal = nDofNod*(mesh.elements{k}(i1)-1)+i2;
            for j1 = 1:nNodElem
                for j2 = 1:nDofNod
                    iRow = j2+(j1-1)*nDofNod; 
                    iRowGlobal = nDofNod*(mesh.elements{k}(j1)-1)+j2;
                    
                    iloc = iloc+1;
                    if iloc > nSparse
                        % Aumentar tamaño si la estimación fue corta
                        nSparse = nSparse * 2;
                        iRowSparse(nSparse, 1) = 0;
                        iColSparse(nSparse, 1) = 0;
                        valueSparseK(nSparse, 1) = 0;
                    end
                    iRowSparse(iloc,1) = iRowGlobal; 
                    iColSparse(iloc,1) = iColGlobal; 
                    valueSparseK(iloc,1) = kEle(iRow,iCol);
                end
            end
        end
    end
    
    % Ensamblaje de F
    Fe_mat_body = reshape(Fe_body, nDofNod, nNodElem)';
    Fe_mat_pressure = reshape(Fe_pressure, nDofNod, nNodElem)';
    
    F(mesh.elements{k}, :) = F(mesh.elements{k}, :) + Fe_mat_body + Fe_mat_pressure;
end
% Limpiar exceso de sparse y crear K
iRowSparse = iRowSparse(1:iloc);
iColSparse = iColSparse(1:iloc);
valueSparseK = valueSparseK(1:iloc);
K = sparse(iRowSparse,iColSparse,valueSparseK, nNod*nDofNod, nNod*nDofNod);
% Reshape F
F = reshape(F',[],1);

%% BC
free = true(nNod,nDofNod);
tol = 1e-6; 
xmax = max(mesh.nodes(:,1));
for node_row = 1:nNod % Iterar por FILAS
    % Fijar nodos en la base (y=0)
    if mesh.nodes(node_row, 1) == 0
            free(node_row, 1) = false;
    else
        if abs(mesh.nodes(node_row, 2)) < tol
                free(node_row, 2) = false;
        end
    end
end
free = reshape(free',[],1);

%% Solucion

fprintf('Resolviendo sistema lineal...\n');
D = zeros(size(F,1),1);
D(free) = K(free,free)\F(free);
R = K*D; % Reacciones
% Preparo matrices para postprocesado
D_mat = reshape(D',nDofNod,nNod)';
R_mat = reshape(R',nDofNod,nNod)';
max_disp_x = max(abs(D_mat(:,1)));
max_disp_y = max(abs(D_mat(:,2)));
fprintf('Solución completa.\n');
fprintf('Desplazamiento X máximo: %.4f mm\n', max_disp_x);
fprintf('Desplazamiento Y máximo: %.4f mm\n', max_disp_y);
%% Tensiones

stressDim = size(C,1); % 3 para 2D
npg_stress = npg_vol; % Usar los mismos puntos que para la rigidez
[~, lgp_stress] = gauss1D(npg_stress);

% Chequeo rápido si hay elementos CST
isCST = any(cellfun(@(c) numel(c) == 3, mesh.elements));
isQ4Q8 = any(cellfun(@(c) numel(c) ~= 3, mesh.elements));

if isQ4Q8
    nGP_total = npg_stress^2;
    % S guardará [sigma_x, sigma_y, tau_xy] para cada GP, para cada Elemento
    S = zeros(stressDim, nGP_total, nElem); 
    S_vm = zeros(nGP_total, nElem);
end
if isCST
    S_cst = zeros(stressDim, 1, nElem); % Para CST (tensión constante)
    S_vm_cst = zeros(1, nElem);
end

D_mat_stress = reshape(D,nDofNod,[])'; 

for k = 1:nElem
    elemNodeRows = mesh.elements{k}; 
    localNodes = mesh.nodes(elemNodeRows, :);
    nNodElem = numel(elemNodeRows);
    
    elemDisplacements = reshape(D_mat_stress(elemNodeRows,:)', [], 1);
    
    elemType = '';
    if nNodElem == 3
        elemType = 'CST';
        B = zeros(stressDim, nNodElem * nDofNod);
        
        x1 = localNodes(1,1); y1 = localNodes(1,2);
        x2 = localNodes(2,1); y2 = localNodes(2,2);
        x3 = localNodes(3,1); y3 = localNodes(3,2);
        
        Area = 0.5 * det([1 x1 y1; 1 x2 y2; 1 x3 y3]);
        if Area <= 0
            warning('Elemento %d tiene área 0 o negativa.', k);
            continue;
        end
        
        b1 = y2-y3; b2 = y3-y1; b3 = y1-y2;
        c1 = x3-x2; c2 = x1-x3; c3 = x2-x1;
        
        B(1, 1:2:end) = [b1, b2, b3];
        B(2, 2:2:end) = [c1, c2, c3];
        B(3, 1:2:end) = [c1, c2, c3];
        B(3, 2:2:end) = [b1, b2, b3];
        B = B / (2 * Area);
        
        stress = C * B * elemDisplacements;
        S_cst(:, 1, k) = stress; % Guardar la tensión constante
        
        % Von Mises
        s_x = stress(1); s_y = stress(2); t_xy = stress(3);
        S_vm_cst(1, k) = sqrt(s_x^2 - s_x*s_y + s_y^2 + 3*t_xy^2);
        
    else
        % --- Cálculo para Q4 / Q8 ---
        if nNodElem == 4
            elemType = 'Q4';
        elseif nNodElem == 8
            elemType = 'Q8';
        end

        counterGaussPoint = 0;
        for i = 1:npg_stress
            for j = 1:npg_stress
                counterGaussPoint = counterGaussPoint + 1;
                xi = lgp_stress(i);
                eta = lgp_stress(j);
                
                dN_iso = shapefunsDer([xi eta], elemType);
                Jacobian = dN_iso * localNodes;
                dN_xy = Jacobian \ dN_iso; % dN/dx, dN/dy
                
                B = zeros(stressDim, nNodElem * nDofNod);
                B(1, 1:2:end) = dN_xy(1,:);
                B(2, 2:2:end) = dN_xy(2,:);
                B(3, 1:2:end) = dN_xy(2,:);
                B(3, 2:2:end) = dN_xy(1,:);
                
                stress = C * B * elemDisplacements;
                S(:, counterGaussPoint, k) = stress;
                
                % Von Mises
                s_x = stress(1); s_y = stress(2); t_xy = stress(3);
                S_vm(counterGaussPoint, k) = sqrt(s_x^2 - s_x*s_y + s_y^2 + 3*t_xy^2);
            end
        end
    end
end

if isQ4Q8
    max_vm_q4_q8 = max(S_vm, [], 'all');
    fprintf('Tensión Von Mises Máxima (Q4/Q8): %.4f MPa\n', max_vm_q4_q8);
end
if isCST
    max_vm_cst = max(S_vm_cst, [], 'all');
    fprintf('Tensión Von Mises Máxima (CST): %.4f MPa\n', max_vm_cst);
end
%% Promediar Tensiones para Ploteo 

S_vm_avg = zeros(nElem, 1);
for k = 1:nElem
    nNodElem = numel(mesh.elements{k});
    if nNodElem == 3
        % Para CST, la tensión es constante
        S_vm_avg(k) = S_vm_cst(1, k);
    elseif isQ4Q8
        % Promedio de los puntos de Gauss para Q4/Q8
        S_vm_avg(k) = mean(S_vm(:, k));
    end
end

%% Deformada
scale_factor = 100; 
fprintf('Factor de escala para ploteo: %d\n', scale_factor);

figure;
ax = gca; 
hold(ax, 'on');

plotMalla(ax, mesh, 'b'); % Azul
h_original = findobj(ax, 'Type', 'line', 'Color', 'b'); 

mesh_deformed = mesh;
mesh_deformed.nodes = mesh_deformed.nodes + D_mat * scale_factor;
plotMalla(ax, mesh_deformed, 'r'); % Rojo
h_deformada = findobj(ax, 'Type', 'line', 'Color', 'r'); 

hold(ax, 'off');
title(ax, ['Forma Deformada (Escala: ', num2str(scale_factor), ') - Plane Strain']); 

if ~isempty(h_original) && ~isempty(h_deformada)
    legend(ax, [h_original(1), h_deformada(1)], 'Original', 'Deformada', 'Location', 'northwest');
else
    legend(ax, 'Original', 'Deformada', 'Location', 'northwest'); 
end

axis(ax, 'equal');
xlabel(ax, 'Coordenada X (mm)');
ylabel(ax, 'Coordenada Y (mm)');
grid(ax, 'on');


%% Ploteo de Tensiones 

figure;
ax_stress = gca;
hold(ax_stress, 'on');
title(ax_stress, 'Tensiones de Von Mises (MPa)');
xlabel(ax_stress, 'Coordenada X (mm)');
ylabel(ax_stress, 'Coordenada Y (mm)');
axis(ax_stress, 'equal');
grid(ax_stress, 'on');

max_stress_plot = max(S_vm_avg);
min_stress_plot = min(S_vm_avg);
fprintf('Rango de tensiones (Von Mises): %.4f a %.4f MPa\n', min_stress_plot, max_stress_plot);


for k = 1:nElem
    elemNodeRows = mesh.elements{k};
    localNodes = mesh.nodes(elemNodeRows, :);
    
    plot_nodes_coords = localNodes;
    
    if numel(elemNodeRows) == 8
        plot_nodes_coords = localNodes([1, 2, 3, 4], :);
    end
    
    % Dibujar el 'patch' (parche) del elemento
    % El color se define por S_vm_avg(k)
    patch(ax_stress, plot_nodes_coords(:, 1), plot_nodes_coords(:, 2), S_vm_avg(k), 'EdgeColor', 'none');
end

% Configurar el gráfico
shading(ax_stress, 'flat'); % 'flat' usa un color por elemento
colorbar(ax_stress);
clim(ax_stress, [min_stress_plot, max_stress_plot]); % Escala de color
hold(ax_stress, 'off');

%% Momentos

datos.H_presa   = max(mesh.nodes(:,2))/1000;    % [m]
datos.B_base    = max(mesh.nodes(:,1))/1000;    % [m]
datos.c_cresta  = max(mesh.nodes(mesh.nodes(:,2) == max(mesh.nodes(:,2)),1))/1000 - min(mesh.nodes(mesh.nodes(:,2) == max(mesh.nodes(:,2)),1))/1000;     % [m]
datos.d_libre   = datos.H_presa - water_level_meters;     % [m]
datos.rho_c     = 2400;  % [kg/m^3]
datos.rho_a     = 1000;  % [kg/m^3]

x_pivot_talon_mm = datos.B_base * 1000; % Pivote en [mm]
tol = 1e-6;
base_node_indices = find(abs(mesh.nodes(:, 2)) < tol);

M_reacciones_FEM_Nmm = 0; % En N*mm
for i = 1:length(base_node_indices)
    node_idx = base_node_indices(i);
    Ry = R_mat(node_idx, 2); % Reacción vertical [N]
    x_node = mesh.nodes(node_idx, 1); % Posición [mm]
    brazo_palanca = x_pivot_talon_mm - x_node; % [mm]
    M_reacciones_FEM_Nmm = M_reacciones_FEM_Nmm + (Ry * brazo_palanca);
end

M_reacciones_FEM_Nm = M_reacciones_FEM_Nmm / 1000; % Convertir a N*m

[M_estab_Nm, M_desestab_Nm, FS_vuelco] = calcularMomentosVuelco(datos);

M_neto_analitico_Nm = M_estab_Nm - M_desestab_Nm;

error_relativo = abs(abs(M_reacciones_FEM_Nm) - abs(M_neto_analitico_Nm)) / abs(M_neto_analitico_Nm);


fprintf('Momento Reacciones (FEM):      %.3f kNm\n', M_reacciones_FEM_Nm / 1000);
fprintf('Momento Externo Neto (Analítico): %.3f kNm\n', -M_neto_analitico_Nm / 1000);
fprintf('Error Relativo (Equilibrio):   %.4f %%\n', error_relativo * 100);
fprintf('Factor de Seguridad (Vuelco):  %.3f\n', FS_vuelco);

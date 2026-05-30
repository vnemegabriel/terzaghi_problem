function plotMalla(ax, mesh, color)
% plotMalla Dibuja una malla 2D (TRIAs y QUADs) en un eje específico.
%
% VERSIÓN MODIFICADA: Eliminada la dependencia de 'nodeIDMap'.
% Se asume que mesh.elements{i} contiene los ÍNDICES DE FILA
% para mesh.nodes.
%
% SINOPSIS:
%   plotMalla(ax, mesh_struct, 'b')
%
% INPUT:
%   ax   - El handle del axes donde se dibujará (ej: gca)
%   mesh - La estructura generada por 'lectorMalla'.
%   color- String del color (ej: 'b', 'r', 'k-')

% --- Verificación de entrada ---
% Se elimina 'nodeIDMap' de los campos requeridos
requiredFields = {'nodes', 'elements'};
if ~all(isfield(mesh, requiredFields))
    error(['El struct ''mesh'' debe contener los campos: ' ...
           '''nodes'', ''elements'', y ''plane''.']);
end

% --- Configuración del Gráfico ---
% No crea una nueva figura, usa el axes_handle 'ax'
hold(ax, 'on');

% --- Iteración sobre los elementos ---
nElem = length(mesh.elements);
allNodesCoords = mesh.nodes;
% idMap = mesh.nodeIDMap; % <-- ELIMINADO

for i = 1:nElem
    % Se asume que 'mesh.elements{i}' contiene los índices de fila
    rowIndices = mesh.elements{i};
    
    % --- Bloque try/catch eliminado ---
    
    elemCoords = allNodesCoords(rowIndices, :);
    
    % Define el orden de ploteo para cerrar el polígono
    if numel(rowIndices) == 8 % Q8
        % Plotear solo los 4 nodos de esquina para el contorno
        plotOrder = [1, 2, 3, 4, 1];
        plotCoords = elemCoords(plotOrder, :);
    else % Q4 or CST
        plotOrder = [1:numel(rowIndices), 1]; % Cierra el elemento
        plotCoords = elemCoords(plotOrder, :);
    end

    % Dibujar las aristas del elemento en el axes 'ax'
    plot(ax, plotCoords(:, 1), plotCoords(:, 2), 'Color', color, 'LineWidth', 0.5);
end

hold(ax, 'off');


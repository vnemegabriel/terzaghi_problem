function mesh = lectorMalla(filename, plane)
% lectorMalla Lee un archivo .dat de malla (formato Nastran) y extrae nodos y elementos 2D.
%
% VERSIÓN 7.2: Modificada para leer mallas 2D en planos XY, YZ, o XZ.
%              Devuelve el mapa de traducción de IDs de nodo.
%
% SINOPSIS:
%   mesh = lectorMalla('nombre_archivo.dat', 'xy') % Para plano XY
%   mesh = lectorMalla('nombre_archivo.dat', 'yz') % Para plano YZ (predeterminado si se omite)
%   mesh = lectorMalla('nombre_archivo.dat', 'xz') % Para plano XZ
%
% INPUT:
%   filename - Nombre del archivo de malla (ej: 'malla_final.dat')
%   plane    - (Opcional) String especificando el plano 2D: 'xy', 'yz', o 'xz'.
%              Si se omite, el valor predeterminado es 'yz'.
%
% OUTPUT:
%   mesh - Una estructura con los siguientes campos:
%
%     .nodes: [nNod x 2] Matriz con las coordenadas 2D seleccionadas.
%              El orden depende del 'plane' elegido:
%                - 'xy': [X, Y]
%                - 'yz': [Y, Z]
%                - 'xz': [X, Z]
%
%     .elements: {nElem x 1} Array de celdas. Cada celda contiene un
%                vector [N1, N2, N3, (N4)], donde los IDs de nodo
%                son los NÚMEROS ORIGINALES del archivo .dat.
%
%     .nodeIDMap: [containers.Map] Un mapa que traduce el ID original
%                 del nodo (Key) al índice de la fila (Value) donde
%                 se almacena en .nodes. (Ej: map(1001) -> 1)
%
%     .plane: [char] String que indica el plano leído ('xy', 'yz', 'xz').

% --- Procesamiento de Argumentos de Entrada ---
if nargin < 2 || isempty(plane)
    plane = 'yz'; % Valor predeterminado
    fprintf("Advertencia: No se especificó el plano. Se asume 'yz'.\n");
elseif ~ismember(lower(plane), {'xy', 'yz', 'xz'})
    error("El argumento 'plane' debe ser 'xy', 'yz', o 'xz'.");
end
plane = lower(plane); % Asegurar minúsculas

% --- Inicialización ---
fileID = fopen(filename, 'r');
if fileID == -1
    error('Error: No se pudo abrir el archivo: %s', filename);
end
tempNodes = {};         % Almacén temporal para coordenadas
mesh.elements = {};     % Salida 2: Info de elementos
nodeID_to_Index_Map = containers.Map('KeyType', 'double', 'ValueType', 'double');
nodeCounter = 0;
elemCounter = 0;
lineCounter = 0;

% --- Lectura del Archivo ---
while ~feof(fileID)
    line = fgetl(fileID);
    lineCounter = lineCounter + 1;

    if ~ischar(line), break; end

    trimmedLine = strtrim(line);

    % Omitir líneas vacías o comentarios
    if isempty(trimmedLine) || trimmedLine(1) == '$' || ...
       startsWith(trimmedLine, 'BEGIN BULK') || ...
       startsWith(trimmedLine, 'ENDDATA')
        continue;
    end

    % --- Procesamiento de Tarjetas (por ancho fijo) ---

    try
        if startsWith(line, 'GRID*')
            % Formato GRANDE (16-char) con continuación en 2da línea
            if feof(fileID), break; end
            line2 = fgetl(fileID); % Leer la línea de continuación
            lineCounter = lineCounter + 1;

            nodeID = str2double(line(9:24));
            x = str2double(line(41:56)); % <<< LEEMOS X SIEMPRE
            y = str2double(line(57:72));
            z = str2double(line2(9:24));

            if isnan(nodeID) || isnan(x) || isnan(y) || isnan(z) % <<< CHEQUEAMOS LAS 3
                 warning('Línea %d GRID* con NaN, omitida.', lineCounter-1);
                 continue;
            end

            nodeCounter = nodeCounter + 1;

            % Guardar la traducción ID -> índice
            nodeID_to_Index_Map(nodeID) = nodeCounter;

            % --- SELECCIÓN DE COORDENADAS SEGÚN EL PLANO ---
            switch plane
                case 'xy'
                    tempNodes{nodeCounter, 1} = [x, y];
                case 'yz'
                    tempNodes{nodeCounter, 1} = [y, z];
                case 'xz'
                    tempNodes{nodeCounter, 1} = [x, z];
            end

        elseif startsWith(line, 'CQUAD4')
            % Formato PEQUEÑO (8-char).
            currentLine = line;

            while strlength(currentLine) < 56
                if feof(fileID), break; end
                line2 = fgetl(fileID);
                lineCounter = lineCounter + 1;
                currentLine = [currentLine, line2];
            end
            elemID_check = str2double(currentLine(9:16));
            n1 = str2double(currentLine(25:32));
            n2 = str2double(currentLine(33:40));
            n3 = str2double(currentLine(41:48));
            n4 = str2double(currentLine(49:56));

            if isnan(elemID_check) || isnan(n1) || isnan(n2) || isnan(n3) || isnan(n4)
                 warning('Línea %d CQUAD4 con NaN, omitida.', lineCounter-1);
                 continue;
            end

            elemCounter = elemCounter + 1;
            mesh.elements{elemCounter, 1} = [n1, n2, n3, n4];

        elseif startsWith(line, 'CTRIA3')
            % Formato PEQUEÑO (8-char).
            currentLine = line;

            while strlength(currentLine) < 48
                if feof(fileID), break; end
                line2 = fgetl(fileID);
                lineCounter = lineCounter + 1;
                currentLine = [currentLine, line2];
            end

            elemID_check = str2double(currentLine(9:16));
            n1 = str2double(currentLine(25:32));
            n2 = str2double(currentLine(33:40));
            n3 = str2double(currentLine(41:48));

            if isnan(elemID_check) || isnan(n1) || isnan(n2) || isnan(n3)
                 warning('Línea %d CTRIA3 con NaN, omitida.', lineCounter-1);
                 continue;
            end

            elemCounter = elemCounter + 1;
            mesh.elements{elemCounter, 1} = [n1, n2, n3];
        end

    catch ME
        warning('Error procesando línea %d ("%s"). Error: %s. Se omite.', ...
                lineCounter, trimmedLine, ME.message);
    end
end

% --- Finalización ---
fclose(fileID);

% --- Construcción de Salidas Finales ---
% Salida 1: mesh.nodes
if isempty(tempNodes)
    error('No se pudo leer ningún nodo (GRID). Verifica el formato del archivo.');
end
mesh.nodes = cell2mat(tempNodes); % Ahora será [nNod x 2] según el plano

% Añadir el mapa a la salida del struct
mesh.nodeIDMap = nodeID_to_Index_Map;

% Guardar qué plano se leyó
mesh.plane = plane;

% Salida 2: mesh.elements
if elemCounter == 0
    warning('Se leyeron nodos pero 0 elementos (CQUAD4/CTRIA3).');
end

end % Fin de la función lectorMallaend
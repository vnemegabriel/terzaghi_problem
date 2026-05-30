function mesh = generateColumnMesh(W, H, nx, ny, eleType)
% generateColumnMesh  Structured Q4 or Q8 mesh for a W x H column.
%   Nodes numbered left-to-right, bottom-to-top.
%   x in [0, W],  y in [0, H]
%
%   eleType : (optional) 'Q4' (default) or 'Q8'.
%             For Q8, mid-edge nodes are inserted on every element edge.
%
%   mesh.nodes    : [nNod x 2]
%   mesh.elements : {nElem x 1}  each cell [n1 n2 n3 n4 (n5 n6 n7 n8)] CCW

if nargin < 5 || isempty(eleType)
    eleType = 'Q4';
end

validateattributes(W,  {'numeric'}, {'scalar','real','positive','finite'}, mfilename, 'W')
validateattributes(H,  {'numeric'}, {'scalar','real','positive','finite'}, mfilename, 'H')
validateattributes(nx, {'numeric'}, {'scalar','integer','positive'},       mfilename, 'nx')
validateattributes(ny, {'numeric'}, {'scalar','integer','positive'},       mfilename, 'ny')

dx = W / nx;
dy = H / ny;

switch eleType

    % --------------------------------------------------------------- Q4
    case 'Q4'
        nNodX = nx + 1;
        nNodY = ny + 1;
        nNod  = nNodX * nNodY;
        nElem = nx * ny;

        nodes = zeros(nNod, 2);
        for j = 1:nNodY
            for i = 1:nNodX
                idx = (j-1)*nNodX + i;
                nodes(idx, :) = [(i-1)*dx, (j-1)*dy];
            end
        end

        elements = cell(nElem, 1);
        eIdx = 0;
        for j = 1:ny
            for i = 1:nx
                eIdx = eIdx + 1;
                n1 = (j-1)*nNodX + i;
                n2 = (j-1)*nNodX + i + 1;
                n3 =  j   *nNodX + i + 1;
                n4 =  j   *nNodX + i;
                elements{eIdx} = [n1, n2, n3, n4];
            end
        end

    % --------------------------------------------------------------- Q8
    case 'Q8'
        % Two row types alternating from bottom to top:
        %   "full" rows  : corners + horizontal mid-edge nodes (2*nx+1 each)
        %   "vert" rows  : vertical mid-edge nodes only        (nx+1 each)
        nFull = 2*nx + 1;           % nodes per full row
        nVert = nx + 1;             % nodes per vertical-only row
        nNod  = (ny+1)*nFull + ny*nVert;
        nElem = nx * ny;

        nodes = zeros(nNod, 2);

        % Helper offsets for sequential numbering
        % Row block j (j = 0..ny):  full row at y = j*dy,
        % then (if j<ny) a vertical-only row at y = (j+0.5)*dy.
        rowStart = zeros(ny+1, 1);          % start index of each full row
        vertStart = zeros(ny,   1);         % start index of each vertical-only row
        offset = 0;
        for j = 0:ny
            rowStart(j+1) = offset + 1;
            offset = offset + nFull;
            if j < ny
                vertStart(j+1) = offset + 1;
                offset = offset + nVert;
            end
        end

        % Populate full rows (y = j*dy, x = 0, dx/2, dx, 3dx/2, ..., W)
        for j = 0:ny
            y = j*dy;
            for k = 1:nFull
                x = (k-1) * dx/2;
                nodes(rowStart(j+1) + k - 1, :) = [x, y];
            end
        end

        % Populate vertical-only rows (y = (j+0.5)*dy, x = 0, dx, 2dx, ..., W)
        for j = 0:ny-1
            y = (j + 0.5) * dy;
            for k = 1:nVert
                x = (k-1) * dx;
                nodes(vertStart(j+1) + k - 1, :) = [x, y];
            end
        end

        % Element connectivity:
        %   corners (CCW from bottom-left):  n1 n2 n3 n4
        %   mid-edges (bottom, right, top, left): n5 n6 n7 n8
        elements = cell(nElem, 1);
        eIdx = 0;
        for j = 1:ny                       % element row (1-indexed)
            for i = 1:nx                   % element col (1-indexed)
                eIdx = eIdx + 1;

                % Full-row corner indices in row (j-1) and j
                colL = 2*i - 1;            % left  corner col index in full row
                colM = 2*i;                % horizontal mid col index
                colR = 2*i + 1;            % right corner col index

                n1 = rowStart(j)   + colL - 1;     % bottom-left
                n2 = rowStart(j)   + colR - 1;     % bottom-right
                n3 = rowStart(j+1) + colR - 1;     % top-right
                n4 = rowStart(j+1) + colL - 1;     % top-left
                n5 = rowStart(j)   + colM - 1;     % bottom mid
                n7 = rowStart(j+1) + colM - 1;     % top mid

                % Vertical mid-edges in vertical row (j-1) (1-indexed: j)
                % Left mid : vertical col i,   Right mid : vertical col i+1
                n8 = vertStart(j) + i - 1;          % left  mid
                n6 = vertStart(j) + i;              % right mid

                elements{eIdx} = [n1, n2, n3, n4, n5, n6, n7, n8];
            end
        end

    otherwise
        error('generateColumnMesh:invalidEleType', ...
              "eleType must be 'Q4' or 'Q8' (got '%s').", eleType)
end

mesh.nodes    = nodes;
mesh.elements = elements;
mesh.eleType  = eleType;
end

clear 
clc
close all

%Deducted from 6.13 COOK
% 
% mesh.nodes = 1000*[0 0;
%               0.5  0;
%               1 0;
%               0 0.5; 
%               0.5 0.5;
%               1 0.5;
%               0 1;
%               0.5 1;
%               1 1;
%               0 1.5;
%               0.5 1.5;
%               1 1.5];

%Q8 

mesh.nodes = 1000*[0 0;
              0.5  0;
              1 0;
              1.5 0;
              2 0;
              0 0.5; 
              1 0.5;
              2 0.5;
              0 1;
              0.5 1;
              1 1;
              1.5 1;
              2 1;
              0 1.5;
              1 1.5;
              2 1.5;
              0 2;
              0.5 2;
              1 2;
              1.5 2;
              2 2];

% Q4

% mesh.elements = [1 2 5 4;
%                     2 3 6 5;
%                     4 5 8 7;
%                     5 6 9 8;
%                     7 8 11 10;
%                     8 9 12 11];

% Q8

mesh.elements = [1 2 3 7 11 10 9 6;
                 3 4 5 8 13 12 11 7;
                 9 10 11 15 19 18 17 14;
                 11 12 13 16 21 20 19 15];

meshplot(mesh.nodes,mesh.elements)


% Definition of nodes and dofs

nNodElem = 8; nDofNod = 2;
nNod = size(mesh.nodes,1);
nDofTot = nDofNod * nNod;
nElem = size(mesh.elements,1);

%Mats

young = 2e5; pnu=0.33;

%Plane stress constitutive
E = young/(1-pnu^2)*[1 pnu 0;pnu 1 0;0 0 0.5*(1-pnu)];

% Thickness


% Aiming to see a target \sigma_x = 100
sigma_x = 400;

P = 1000; % [N]


F = zeros(nNod,nDofNod);

% sigma x
% Q4
% F(3,1) = P;
% F(6,1) = 2*P;
% F(9,1) = 2*P;
% F(12,1) = P;
% F(4,1) = -2*P;
% F(7,1) = -2*P;
% F(10,1) = -P;

F([5 8 13 16 21],1) = P * [1 1 2 1 1];
F([6 9 14 17],1) = -P * [1 2 1 1];

% Agrego caso sigma_y

% F(3,2) = P;
% F(10,2) = -P;
% F(11,2) = -2*P;
% F(12,2) = -P;

% Corte (Cauchy)
% 
% F(2,1)=-2*P;
% F(3,[1 2])=[-P P];
% F(4,2)=-2*P;
% F(7,2)=-2*P;
% F(6,2)=2*P;
% F(9,2)=2*P;
% F(10,[1 2])=[P -P];
% F(11,1)=2*P;
% F(12, [1 2])=[P P];
% 
% F = reshape(F',[],1);

% BC s

free = true(nNod,nDofNod);

free(1,:) = [0 0];
free(2,2) = 0;

free = reshape(free',[],1);

iloc=0;

for k = 1:nElem
    localNodes = mesh.nodes(mesh.elements(k,:),:);
    Kelem = getStiffnessMatrix(localNodes,'Q8',1,E,[2 2]);
    %%% Assembly
    for i1=1:nNodElem
        for i2=1:nDofNod
            iCol = i2+(i1-1)*nDofNod;
            iColGlobal = nDofNod*(mesh.elements(k,i1)-1)+i2;
            for j1=1:nNodElem
                for j2=1:nDofNod
                    iRow = j2+(j1-1)*nDofNod; 
                    iRowGlobal = nDofNod*(mesh.elements(k,j1)-1)+j2;
                    iloc = iloc+1;
                    iRowSparse(iloc,1) = iRowGlobal; iColSparse(iloc,1) = iColGlobal; 
                    valueSparseK(iloc,1) = Kelem(iRow,iCol);
                end
            end
        end
    end
end

KS = sparse(iRowSparse,iColSparse,valueSparseK); 

KS = KS(free,free);

F = F(free);

Dr=KS\F;

D=zeros(nDofNod*nNod,1);

D(free)=D(free)+Dr;

Dn=reshape(D,nDofNod,nNod)';

S = stressesAtPoints(full(D),mesh,'Q8',E,[2 2]);




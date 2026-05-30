function [M_estab, M_desestab, FS_vuelco] = calcularMomentosVuelco(datos)
% Calcula momentos y FS de vuelco 
% Entradas: struct 'datos' con H_presa, B_base, c_cresta, d_libre, rho_c, rho_a (en metros y kg/m^3)
% Salidas: Momentos en N*m y FS.

g = 9.81; 
t = 0.001;    % Espesor unitario (1mm)

% 1. Momento Desestabilizante (N*m)
H_agua = datos.H_presa - datos.d_libre;
M_desestab = (1/6) * datos.rho_a * g * H_agua^3 * t;

% 2. Momento Estabilizante (N*m)
% Pivote en (x = B_base, y = 0)

% Rectángulo
Area_rect = datos.c_cresta * datos.H_presa;
Peso_rect = Area_rect * t * datos.rho_c * g;
x_rect = datos.c_cresta / 2;
d_palanca_rect = datos.B_base - x_rect;
M_rect = Peso_rect * d_palanca_rect;

% Triángulo
base_tri = datos.B_base - datos.c_cresta;
Area_tri = 0.5 * base_tri * datos.H_presa;
Peso_tri = Area_tri * t * datos.rho_c * g;
x_tri = datos.c_cresta + (1/3) * base_tri;
d_palanca_tri = datos.B_base - x_tri;
M_tri = Peso_tri * d_palanca_tri;

M_estab = M_rect + M_tri;

% 3. Factor de Seguridad
FS_vuelco = M_estab / M_desestab;

end
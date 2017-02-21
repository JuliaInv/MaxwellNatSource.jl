

% syms  Ex1 Ex2  Ey1 Ey2 
% syms  Hx1 Hx2  Hy1 Hy2
% E = [ Ex1 Ex2 ; Ey1 Ey2 ];
% H = [ Hx1 Hx2 ; Hy1 Hy2 ];
% ZZ = E / H ;


syms Ex1r Ex2r Ey1r Ey2r  real
syms Ex1i Ex2i Ey1i Ey2i  real

syms Hx1r Hx2r Hy1r Hy2r  real
syms Hx1i Hx2i Hy1i Hy2i  real

syms Z11 Z21 Z12 Z22

Ex1 = Ex1r + 1i*Ex1i ;
Ex2 = Ex2r + 1i*Ex2i ;

Ey1 = Ey1r + 1i*Ey1i ;
Ey2 = Ey2r + 1i*Ey2i ;

Hx1 = Hx1r + 1i*Hx1i ;
Hx2 = Hx2r + 1i*Hx2i ;

Hy1 = Hy1r + 1i*Hy1i ;
Hy2 = Hy2r + 1i*Hy2i ;

HH = Hx1*Hy2 - Hx2*Hy1;

Z11 = (  Ex1*Hy2 - Ex2*Hy1 ) / HH;
Z21 = (  Ey1*Hy2 - Ey2*Hy1 ) / HH;
Z12 = ( -Ex1*Hx2 + Ex2*Hx1 ) / HH;
Z22 = ( -Ey1*Hx2 + Ey2*Hx1 ) / HH;


Z11r = simplify(expand(real(Z11))) ;
Z11i = simplify(expand(imag(Z11))) ;

Z21r = simplify(expand(real(Z21))) ;
Z21i = simplify(expand(imag(Z21))) ;

Z12r = simplify(expand(real(Z12))) ;
Z12i = simplify(expand(imag(Z12))) ;

Z22r = simplify(expand(real(Z22))) ;
Z22i = simplify(expand(imag(Z22))) ;


ZZ = [ Z11r Z11i  Z21r Z21i  Z12r Z12i  Z22r Z22i ];
fields = [ Ex1r Ex1i  Ey1r Ey1i  Hx1r Hx1i  Hy1r Hy1i  ...
           Ex2r Ex2i  Ey2r Ey2i  Hx2r Hx2i  Hy2r Hy2i ];


ndata   = length(ZZ);
nfields = length(fields);

DD = sym('DD', [ndata  nfields] );
for iz = 1:ndata        
   for ifld = 1:nfields
      DD(iz,ifld) = diff(ZZ(iz), fields(ifld)) ;
   end
end
        
DD = simplify(DD);

% Save to a file:
% matlabFunction(DD,'file', 'MTderivs.m')

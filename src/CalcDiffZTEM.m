

% syms  Hx1 Hx2  Hy1 Hy2   Hz1 Hz2
% H = [ Hx1 Hy1 ; Hx2 Hy2 ];
% Hz = [ Hz1 ; Hz2 ];
% T = H \ Hz ;



syms Hx1r Hx2r Hy1r Hy2r Hz1r Hz2r real
syms Hx1i Hx2i Hy1i Hy2i Hz1i Hz2i real

syms Tx Ty

Hx1 = Hx1r + 1i*Hx1i ;
Hx2 = Hx2r + 1i*Hx2i ;

Hy1 = Hy1r + 1i*Hy1i ;
Hy2 = Hy2r + 1i*Hy2i ;

Hz1 = Hz1r + 1i*Hz1i ;
Hz2 = Hz2r + 1i*Hz2i ;


HH = Hx1*Hy2 - Hx2*Hy1;

Tx = ( -Hy1*Hz2 + Hy2*Hz1 ) / HH;
Ty = (  Hx1*Hz2 - Hx2*Hz1 ) / HH;


Txr = simplify(expand(real(Tx))) ;
Txi = simplify(expand(imag(Tx))) ;

Tyr = simplify(expand(real(Ty))) ;
Tyi = simplify(expand(imag(Ty))) ;


TT = [ Txr Txi  Tyr Tyi ];
fields = [  Hx1r Hx1i  Hy1r Hy1i  Hz1r Hz1i  ...
            Hx2r Hx2i  Hy2r Hy2i  Hz2r Hz2i];


ndata   = length(TT);
nfields = length(fields);

DD = sym('DD', [ndata  nfields] );
for iz = 1:ndata        
   for ifld = 1:nfields
      DD(iz,ifld) = diff(TT(iz), fields(ifld)) ;
   end
end
 
DD = simplify(DD);

% Save to a file:
% matlabFunction(DD,'file', 'ZTEMderivs.m')

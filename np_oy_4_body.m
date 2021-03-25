x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];

J = 1;
x_113  = sparse(kron(kron(eye(2^0),x),eye(2^13)));
z_113  = sparse(kron(kron(eye(2^0),z),eye(2^13)));

x_122  = sparse(kron(kron(eye(2^1),x),eye(2^12)));
z_122  = sparse(kron(kron(eye(2^1),z),eye(2^12)));

x_114  = sparse(kron(kron(eye(2^2),x),eye(2^11)));
z_114  = sparse(kron(kron(eye(2^2),z),eye(2^11)));

x_124  = sparse(kron(kron(eye(2^3),x),eye(2^10)));

x_211  = sparse(kron(kron(eye(2^4),x),eye(2^9)));
z_211  = sparse(kron(kron(eye(2^4),z),eye(2^9)));

x_221  = sparse(kron(kron(eye(2^5),x),eye(2^8)));

x_213  = sparse(kron(kron(eye(2^6),x),eye(2^7)));

x_222  = sparse(kron(kron(eye(2^7),x),eye(2^6)));

z_111  = sparse(kron(kron(eye(2^8),z),eye(2^5)));

z_112  = sparse(kron(kron(eye(2^9),z),eye(2^4)));

z_014  = sparse(kron(kron(eye(2^10),z),eye(2^3)));

z_103  = sparse(kron(kron(eye(2^11),z),eye(2^2)));

x_a111 = sparse(kron(kron(eye(2^12),x),eye(2^1)));
y_a111 = sparse(kron(kron(eye(2^12),y),eye(2^1)));
z_a111 = sparse(kron(kron(eye(2^12),z),eye(2^1)));

x_a112 = sparse(kron(kron(eye(2^13),x),eye(2^0)));
y_a112 = sparse(kron(kron(eye(2^13),y),eye(2^0)));
z_a112 = sparse(kron(kron(eye(2^13),z),eye(2^0)));

I = sparse(eye(2^14));

A = 1; U = 1; t = 1;

H_8_body = sparse(-J*(x_113*x_122*x_114*x_124*x_211*x_221*x_213*x_222 + z_111*z_112*z_113*z_114 + z_014*z_111 + z_112*z_103 + z_114*z_211 + z_113*z_122));

H_4_body = sparse(-A*(z_111*z_112*z_113*z_114 + z_014*z_111 + z_112*z_103 + z_114*z_211 + z_113*z_122 ...
    + (I - z_a111 + z_a112 + z_a111*z_a112) * (z_a111 + z_a112 + z_a111*z_a112 - I) ...
    + (I + z_a111 - z_a112 + z_a111*z_a112) * (I - z_a111 - z_a112 - z_a111*z_a112)) ...
    + (U/2)*(z_a111 + z_a112 + z_a111*z_a112 - I) ...
    + (t/2)*((x_a112 + z_a111*x_a112)*x_113*x_114 + (x_a111*x_a112 + y_a111*y_a112)*x_122*x_124 ...
    + (x_a112 - z_a111*x_a112)*x_221*x_222 + (x_a111*x_a112 - y_a111*y_a112)*x_211*x_213 ));

abs(eigs(H_8_body, 1, 'smallestreal')-eigs(H_4_body, 1, 'smallestreal'))

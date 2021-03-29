x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];

J = 1;
x_113  = sparse(kron(kron(speye(2^0),x),speye(2^17)));
z_113  = sparse(kron(kron(speye(2^0),z),speye(2^17)));

x_122  = sparse(kron(kron(speye(2^1),x),speye(2^16)));
z_122  = sparse(kron(kron(speye(2^1),z),speye(2^16)));

x_114  = sparse(kron(kron(speye(2^2),x),speye(2^15)));
z_114  = sparse(kron(kron(speye(2^2),z),speye(2^15)));

x_124  = sparse(kron(kron(speye(2^3),x),speye(2^14)));

x_211  = sparse(kron(kron(speye(2^4),x),speye(2^13)));
z_211  = sparse(kron(kron(speye(2^4),z),speye(2^13)));

x_221  = sparse(kron(kron(speye(2^5),x),speye(2^12)));

x_213  = sparse(kron(kron(speye(2^6),x),speye(2^11)));

x_222  = sparse(kron(kron(speye(2^7),x),speye(2^10)));

z_111  = sparse(kron(kron(speye(2^8),z),speye(2^9)));

z_112  = sparse(kron(kron(speye(2^9),z),speye(2^8)));

z_014  = sparse(kron(kron(speye(2^10),z),speye(2^7)));

z_103  = sparse(kron(kron(speye(2^11),z),speye(2^6)));

x_a111 = sparse(kron(kron(speye(2^12),x),speye(2^5)));
y_a111 = sparse(kron(kron(speye(2^12),y),speye(2^5)));
z_a111 = sparse(kron(kron(speye(2^12),z),speye(2^5)));

x_a112 = sparse(kron(kron(speye(2^13),x),speye(2^4)));
y_a112 = sparse(kron(kron(speye(2^13),y),speye(2^4)));
z_a112 = sparse(kron(kron(speye(2^13),z),speye(2^4)));

x_a211 = sparse(kron(kron(speye(2^14),x),speye(2^3)));
y_a211 = sparse(kron(kron(speye(2^14),y),speye(2^3)));
z_a211 = sparse(kron(kron(speye(2^14),z),speye(2^3)));

x_a212 = sparse(kron(kron(speye(2^15),x),speye(2^2)));
y_a212 = sparse(kron(kron(speye(2^15),y),speye(2^2)));
z_a212 = sparse(kron(kron(speye(2^15),z),speye(2^2)));

x_a121 = sparse(kron(kron(speye(2^16),x),speye(2^1)));
y_a121 = sparse(kron(kron(speye(2^16),y),speye(2^1)));
z_a121 = sparse(kron(kron(speye(2^16),z),speye(2^1)));

x_a122 = sparse(kron(kron(speye(2^17),x),speye(2^0)));
y_a122 = sparse(kron(kron(speye(2^17),y),speye(2^0)));
z_a122 = sparse(kron(kron(speye(2^17),z),speye(2^0)));

I = speye(2^18);

U = 1; t = 1;

H_8_body = sparse(-J*(x_113*x_122*x_114*x_124*x_211*x_221*x_213*x_222 + z_111*z_112*z_113*z_114 + z_014*z_111 + z_112*z_103 + z_114*z_211 + z_113*z_122));

H_4_body = sparse(J*(-z_111*z_112*z_113*z_114 - z_014*z_111 + z_112*z_103 + z_114*z_211 + z_113*z_122 ...
    + (I - z_a111 + z_a112 + z_a111*z_a112) * (z_a121 + z_a122 + z_a121*z_a122 - I) ...
    + (I + z_a111 - z_a112 + z_a111*z_a112) * (I - z_a211 - z_a212 - z_a211*z_a212)) ...
    - (U/2)*(z_a111 + z_a112 + z_a111*z_a112 - I) ...
    - (t/2)*((x_a112 + z_a111*x_a112)*x_113*x_114 + (x_a111*x_a112 + y_a111*y_a112)*x_122*x_124 ...
    + (x_a112 - z_a111*x_a112)*x_221*x_222 + (x_a111*x_a112 - y_a111*y_a112)*x_211*x_213 ));

abs(eigs(H_8_body, 1, 'smallestreal')-eigs(H_4_body, 1, 'smallestreal'))

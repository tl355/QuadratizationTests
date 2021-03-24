x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];

J = 1;
x_113  = kron(kron(eye(2^0),x),eye(2^17));
z_113  = kron(kron(eye(2^0),z),eye(2^17));

x_122  = kron(kron(eye(2^1),x),eye(2^16));
z_122  = kron(kron(eye(2^1),z),eye(2^16));

x_114  = kron(kron(eye(2^2),x),eye(2^15));
z_114  = kron(kron(eye(2^2),z),eye(2^15));

x_124  = kron(kron(eye(2^3),x),eye(2^14));

x_211  = kron(kron(eye(2^4),x),eye(2^13));
z_211  = kron(kron(eye(2^4),z),eye(2^13));

x_221  = kron(kron(eye(2^5),x),eye(2^12));

x_213  = kron(kron(eye(2^6),x),eye(2^11));

x_222  = kron(kron(eye(2^7),x),eye(2^10));

z_111  = kron(kron(eye(2^8),z),eye(2^9));

z_112  = kron(kron(eye(2^9),z),eye(2^8));

z_014  = kron(kron(eye(2^10),z),eye(2^7));

z_103  = kron(kron(eye(2^11),z),eye(2^6));


x_a111 = kron(kron(eye(2^12),x),eye(2^5));
y_a111 = kron(kron(eye(2^12),y),eye(2^5));
z_a111 = kron(kron(eye(2^12),z),eye(2^5));

x_a112 = kron(kron(eye(2^13),x),eye(2^4));
y_a112 = kron(kron(eye(2^13),y),eye(2^4));
z_a112 = kron(kron(eye(2^13),z),eye(2^4));

z_a211 = kron(kron(eye(2^14),z),eye(2^3));

z_a212 = kron(kron(eye(2^15),z),eye(2^2));

z_a121 = kron(kron(eye(2^16),z),eye(2^1));

z_a122 = kron(kron(eye(2^17),z),eye(2^0));

I_size = 2^18;

A = 1; U = 1; t = 1;

H_8_body = -J*(x_113*x_122*x_114*x_124*x_211*x_221*x_213*x_222 + z_111*z_112*z_113*z_114 + z_014*z_111 + z_112*z_103 + z_114*z_211 + z_113*z_122);

H_4_body = -A*(z_111*z_112*z_113*z_114 + z_014*z_111 + z_112*z_103 + z_114*z_211 + z_113*z_122 ...
    + (eye(I_size) - z_a111 + z_a112 + z_a111*z_a112) * (z_a121 + z_a122 + z_a121*z_a122 - eye(I_size)) ...
    + (eye(I_size) + z_a111 - z_a112 + z_a111*z_a112) * (eye(I_size) - z_a211 - z_a212 - z_a211*z_a212)) ...
    + (U/2)*(z_a111 + z_a112 + z_a111*z_a112 - eye(I_size)) ...
    + (t/2)*((x_a112 + z_a111*x_a112)*x_113*x_114 + (x_a111*x_a112 + y_a111*y_a112)*x_122*x_124 ...
    + (x_a112 - z_a111*x_a112)*x_221*x_222 + (x_a111*x_a112 - y_a111*y_a112)*x_211*x_213 );

abs(min(eig(H_8_body))-min(eig(H_4_body)))
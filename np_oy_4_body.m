x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];

J = 1;
x_113  = kron(kron(eye(1),x),eye(8192));
z_113  = kron(kron(eye(1),z),eye(8192));

x_122  = kron(kron(eye(2),x),eye(4096));
z_122  = kron(kron(eye(2),z),eye(4096));

x_114  = kron(kron(eye(4),x),eye(2048));
z_114  = kron(kron(eye(4),z),eye(2048));

x_124  = kron(kron(eye(8),x),eye(1024));

x_211  = kron(kron(eye(16),x),eye(512));
z_211  = kron(kron(eye(16),z),eye(512));

x_221  = kron(kron(eye(32),x),eye(256));

x_213  = kron(kron(eye(64),x),eye(128));

x_222  = kron(kron(eye(128),x),eye(64));

z_111  = kron(kron(eye(256),z),eye(32));
z_112  = kron(kron(eye(512),z),eye(16));

z_014  = kron(kron(eye(1024),z),eye(8));

z_103  = kron(kron(eye(2048),z),eye(4));


x_a111 = kron(kron(eye(4096),x),eye(2));
y_a111 = kron(kron(eye(4096),y),eye(2));
z_a111 = kron(kron(eye(4096),z),eye(2));
x_a112 = kron(kron(eye(8192),x),eye(1));
y_a112 = kron(kron(eye(8192),y),eye(1));
z_a112 = kron(kron(eye(8192),z),eye(1));

I_size = 2^14;

A = 1; U = 1; t = 1;

H_8_body = -J*(x_113*x_122*x_114*x_124*x_211*x_221*x_213*x_222 + z_111*z_112*z_113*z_114 + z_014*z_111 + z_112*z_103 + z_114*z_211 + z_113*z_122);
H_4_body = -A*(z_111*z_112*z_113*z_114 + z_014*z_111 + z_112*z_103 + z_114*z_211 + z_113*z_122 ...
    + (eye(I_size) - z_a111 + z_a112 + z_a111*z_a112) * (z_a111 + z_a112 + z_a111*z_a112 - eye(I_size)) ...
    + (eye(I_size) + z_a111 - z_a112 + z_a111*z_a112) * (eye(I_size) - z_a111 - z_a112 - z_a111*z_a112)) ...
    + (U/2)*(z_a111 + z_a112 + z_a111*z_a112 - eye(I_size)) ...
    + (t/2)*((x_a112 + z_a111*x_a112)*x_113*x_114 + (x_a111*x_a112 + y_a111*y_a112)*x_122*x_124 ...
    + (x_a112 - z_a111*x_a112)*x_221*x_222 + (x_a111*x_a112 - y_a111*y_a112)*x_211*x_213 );

abs(min(eig(H_8_body))-min(eig(H_4_body)));
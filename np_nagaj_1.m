x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];
gm1 = [0 1 0 ; 1 0 0 ; 0 0 0]; gm2 = [0 -1i 0 ; 1i 0 0 ; 0 0 0]; gm3 = [1 0 0 ; 0 -1 0 ; 0 0 0]; 
gm4 = [0 0 1 ; 0 0 0 ; 1 0 0]; gm5 = [0 0 -1i ; 0 0 0 ; 1i 0 0]; gm6 = [0 0 0 ; 0 0 1 ; 0 1 0];
gm7 = [0 0 0 ; 0 0 -1i ; 0 1i 0]; gm8 = (1/sqrt(3))*[1 0 0 ; 0 1 0 ; 0 0 -2];

x1 = kron(kron(speye(1),x),speye(128));
y1 = kron(kron(speye(1),y),speye(128));
z1 = kron(kron(speye(1),z),speye(128));

x2 = kron(kron(speye(2),x),speye(64));
y2 = kron(kron(speye(2),y),speye(64));
z2 = kron(kron(speye(2),z),speye(64));

x3 = kron(kron(speye(4),x),speye(32));
y3 = kron(kron(speye(4),y),speye(32));
z3 = kron(kron(speye(4),z),speye(32));

x4 = kron(kron(speye(8),x),speye(16));
y4 = kron(kron(speye(8),y),speye(16));
z4 = kron(kron(speye(8),z),speye(16));

xa1 = kron(kron(speye(16),x),speye(8));
ya1 = kron(kron(speye(16),y),speye(8));
za1 = kron(kron(speye(16),z),speye(8));

xa2 = kron(kron(speye(32),x),speye(4));
ya2 = kron(kron(speye(32),y),speye(4));
za2 = kron(kron(speye(32),z),speye(4));

xa3 = kron(kron(speye(64),x),speye(2));
ya3 = kron(kron(speye(64),y),speye(2));
za3 = kron(kron(speye(64),z),speye(2));

xa4 = kron(kron(speye(128),x),speye(1));
ya4 = kron(kron(speye(128),y),speye(1));
za4 = kron(kron(speye(128),z),speye(1));

lam_1a5 = kron(kron(speye(256),gm1),speye(1));
lam_2a5 = kron(kron(speye(256),gm2),speye(1));
lam_4a5 = kron(kron(speye(256),gm4),speye(1));
lam_5a5 = kron(kron(speye(256),gm5),speye(1));
lam_6a5 = kron(kron(speye(256),gm6),speye(1));

lam_1a6 = kron(kron(speye(512),gm1),speye(1));
lam_2a6 = kron(kron(speye(512),gm2),speye(1));
lam_4a6 = kron(kron(speye(512),gm4),speye(1));
lam_5a6 = kron(kron(speye(512),gm5),speye(1));
lam_6a6 = kron(kron(speye(512),gm6),speye(1));

lam_1a7 = kron(kron(speye(1024),gm1),speye(1));
lam_2a7 = kron(kron(speye(1024),gm2),speye(1));
lam_4a7 = kron(kron(speye(1024),gm4),speye(1));
lam_5a7 = kron(kron(speye(1024),gm5),speye(1));
lam_6a7 = kron(kron(speye(1024),gm6),speye(1));

lam_1a8 = kron(kron(speye(2048),gm1),speye(1));
lam_2a8 = kron(kron(speye(2048),gm2),speye(1));
lam_4a8 = kron(kron(speye(2048),gm4),speye(1));
lam_5a8 = kron(kron(speye(2048),gm5),speye(1));
lam_6a8 = kron(kron(speye(2048),gm6),speye(1));

lam_4a9 = kron(kron(speye(4096),gm4),speye(1));
lam_5a9 = kron(kron(speye(4096),gm5),speye(1));
lam_6a9 = kron(kron(speye(4096),gm6),speye(1));

lam_1a10 = kron(kron(speye(8192),gm1),speye(1));
lam_2a10 = kron(kron(speye(8192),gm2),speye(1));
lam_4a10 = kron(kron(speye(8192),gm4),speye(1));
lam_5a10 = kron(kron(speye(8192),gm5),speye(1));
lam_6a10 = kron(kron(speye(8192),gm6),speye(1));

LHS = (1/4)*(x1*x2 + y1*y2 + x1*x2*z3 + x1*x2*x4 + y1*y2*z3 + y1*y2*x4 - x1*x2*z3*x4 - y1*y2*z3*x4);
RHS = (1/2)*(2*lam_6a8 + x1*lam_1a5 + y1*lam_2a5 + lam_6a5 - z3*lam_6a5 + xa1*lam_4a5 + ya1*lam_5a5 + xa1*lam_1a6 ...
    + ya1*lam_2a6 + 2*x4*lam_6a6 + xa2*lam_4a6 + ya2*lam_5a6 + xa2*lam_1a7 + ya2*lam_2a7 + lam_6a7 - za1*lam_6a7 ...
    + x2*lam_4a7 + y2*lam_5a7 + x1*lam_1a8 + y1*lam_2a8 + z3*lam_6a8 + xa3*lam_4a9 + ya3*lam_5a9 + lam_6a9 ...
    + xa4*lam_4a9 + ya4*lam_5a9 + xa4*lam_1a10 + ya4*lam_2a10 + lam_6a10 + z3*lam_6a10 + x2*lam_4a10 + y2*lam_5a10);

abs(min(eig(LHS))-min(eig(RHS)))

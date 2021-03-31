x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];
gm1 = [0 1 0 ; 1 0 0 ; 0 0 0]; gm2 = [0 -1i 0 ; 1i 0 0 ; 0 0 0]; gm3 = [1 0 0 ; 0 -1 0 ; 0 0 0]; 
gm4 = [0 0 1 ; 0 0 0 ; 1 0 0]; gm5 = [0 0 -1i ; 0 0 0 ; 1i 0 0]; gm6 = [0 0 0 ; 0 0 1 ; 0 1 0];
gm7 = [0 0 0 ; 0 0 -1i ; 0 1i 0]; gm8 = (1/sqrt(3))*[1 0 0 ; 0 1 0 ; 0 0 -2];

x1 = kron(kron(speye(2^0),x),speye(2^7));
y1 = kron(kron(speye(2^0),y),speye(2^7));
z1 = kron(kron(speye(2^0),z),speye(2^7));

x2 = kron(kron(speye(2^1),x),speye(2^6));
y2 = kron(kron(speye(2^1),y),speye(2^6));
z2 = kron(kron(speye(2^1),z),speye(2^6));

x3 = kron(kron(speye(2^2),x),speye(2^5));
y3 = kron(kron(speye(2^2),y),speye(2^5));
z3 = kron(kron(speye(2^2),z),speye(2^5));

x4 = kron(kron(speye(2^3),x),speye(2^4));
y4 = kron(kron(speye(2^3),y),speye(2^4));
z4 = kron(kron(speye(2^3),z),speye(2^4));

xa1 = kron(kron(speye(2^4),x),speye(2^3));
ya1 = kron(kron(speye(2^4),y),speye(2^3));
za1 = kron(kron(speye(2^4),z),speye(2^3));

xa2 = kron(kron(speye(2^5),x),speye(2^2));
ya2 = kron(kron(speye(2^5),y),speye(2^2));
za2 = kron(kron(speye(2^5),z),speye(2^2));

xa5 = kron(kron(speye(2^6),x),speye(2^1));
ya5 = kron(kron(speye(2^6),y),speye(2^1));
za5 = kron(kron(speye(2^6),z),speye(2^1));

xa6 = kron(kron(speye(2^7),x),speye(2^0));
ya6 = kron(kron(speye(2^7),y),speye(2^0));
za6 = kron(kron(speye(2^7),z),speye(2^0));

lam_1a5 = kron(kron(speye(3^0),gm1),speye(3^5));
lam_2a5 = kron(kron(speye(3^0),gm2),speye(3^5));
lam_4a5 = kron(kron(speye(3^0),gm4),speye(3^5));
lam_5a5 = kron(kron(speye(3^0),gm5),speye(3^5));
lam_6a5 = kron(kron(speye(3^0),gm6),speye(3^5));

lam_1a6 = kron(kron(speye(3^1),gm1),speye(3^4));
lam_2a6 = kron(kron(speye(3^1),gm2),speye(3^4));
lam_4a6 = kron(kron(speye(3^1),gm4),speye(3^4));
lam_5a6 = kron(kron(speye(3^1),gm5),speye(3^4));
lam_6a6 = kron(kron(speye(3^1),gm6),speye(3^4));

lam_1a7 = kron(kron(speye(3^2),gm1),speye(3^3));
lam_2a7 = kron(kron(speye(3^2),gm2),speye(3^3));
lam_4a7 = kron(kron(speye(3^2),gm4),speye(3^3));
lam_5a7 = kron(kron(speye(3^2),gm5),speye(3^3));
lam_6a7 = kron(kron(speye(3^2),gm6),speye(3^3));

lam_1a8 = kron(kron(speye(3^3),gm1),speye(3^2));
lam_2a8 = kron(kron(speye(3^3),gm2),speye(3^2));
lam_4a8 = kron(kron(speye(3^3),gm4),speye(3^2));
lam_5a8 = kron(kron(speye(3^3),gm5),speye(3^2));
lam_6a8 = kron(kron(speye(3^3),gm6),speye(3^2));

lam_4a9 = kron(kron(speye(3^4),gm4),speye(3^1));
lam_5a9 = kron(kron(speye(3^4),gm5),speye(3^1));
lam_6a9 = kron(kron(speye(3^4),gm6),speye(3^1));

lam_1a10 = kron(kron(speye(3^5),gm1),speye(3^0));
lam_2a10 = kron(kron(speye(3^5),gm2),speye(3^0));
lam_4a10 = kron(kron(speye(3^5),gm4),speye(3^0));
lam_5a10 = kron(kron(speye(3^5),gm5),speye(3^0));
lam_6a10 = kron(kron(speye(3^5),gm6),speye(3^0));

LHS = (1/4)*(kron(x1*x2,eye(3^6)) + kron(y1*y2,eye(3^6)) + kron(x1*x2*z3,eye(3^6)) + kron(x1*x2*x4,eye(3^6)) + kron(y1*y2*z3,eye(3^6)) + kron(y1*y2*x4,eye(3^6)) - kron(x1*x2*z3*x4,eye(3^6)) - kron(y1*y2*z3*x4,eye(3^6)));
RHS = (1/2)*(2*kron(eye(2^8),lam_6a8) + kron(x1,lam_1a5) + kron(y1,lam_2a5) + kron(eye(2^8),lam_6a5) - kron(z3,lam_6a5) + kron(xa1,lam_4a5) + kron(ya1,lam_5a5) + kron(xa1,lam_1a6) ...
    + kron(ya1,lam_2a6) + 2*kron(x4,lam_6a6) + kron(xa2,lam_4a6) + kron(ya2,lam_5a6) + kron(xa2,lam_1a7) + kron(ya2,lam_2a7) + kron(eye(2^8),lam_6a7) - kron(za1,lam_6a7) ...
    + kron(x2,lam_4a7) + kron(y2,lam_5a7) + kron(x1,lam_1a8) + kron(y1,lam_2a8) + kron(z3,lam_6a8) + kron(xa5,lam_4a9) + kron(ya5,lam_5a9) + kron(eye(2^8),lam_6a9) ...
    + kron(xa6,lam_4a9) + kron(ya6,lam_5a9) + kron(xa6,lam_1a10) + kron(ya6,lam_2a10) + kron(eye(2^8),lam_6a10) + kron(z3,lam_6a10) + kron(x2,lam_4a10) + kron(y2,lam_5a10));

abs(eigs(LHS, 1, 'smallestreal')-eigs(RHS, 1, 'smallestreal')) % 7.9568


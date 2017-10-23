syms('v','x','y','u','v','real');

U = [(x^2-1)*(y^2-1);(x^2-1)*(y^2-1)];

Uxx = diff(diff(U(1),x),x);
Uxy = diff(diff(U(1),x),y);
Uyy = diff(diff(U(1),y),y);
Vxx = diff(diff(U(2),x),x);
Vxy = diff(diff(U(2),x),y);
Vyy = diff(diff(U(2),y),y);

gradS = [Uxx + (v+1)/2*Vxy + (1-v)/2*Uyy;
        Vyy + (v+1)/2*Uxy + (1-v)/2*Vxx];
gradS = simplify(gradS);
disp(expand(gradS))






%  multiplies complex vector by a phase factor
%  to minimize rms imaginary parts of selected components
%  routine determines phase factor automatically ...
%  Usage: u = chngph(u,ind);

function u = chngph(u,ind)

x = [ real(u(ind)), imag(u(ind)) ];
S = x'*x;
[U,D] = eig(S);
[D,I] = sort(diag(D));
U1 = U(:,I(2));
v = U1(1)-i*U1(2);
% Note that there is a sign indeterminacy!!
u = v*u;
return
end 

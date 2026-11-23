function psi = legendre_basis(z, count)
%LEGENDRE_BASIS Orthonormal Legendre basis for U[-1,1].
z = z(:); psi = zeros(numel(z), count); psi(:,1) = 1;
if count < 2, return; end
previous = ones(size(z)); current = z; psi(:,2) = sqrt(3)*current;
for degree = 2:count-1
    next = ((2*degree-1)*z.*current-(degree-1)*previous)/degree;
    psi(:,degree+1) = sqrt(2*degree+1)*next;
    previous = current; current = next;
end
end

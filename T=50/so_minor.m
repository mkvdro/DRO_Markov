function val = so_minor(A,d,j,l,k)


M1=A([1:(j-1) (j+1):(d-1)],[1:(l-1) (l+1):d]);
val=M1(:,[1:(k-2) k:(d-1)]);
end


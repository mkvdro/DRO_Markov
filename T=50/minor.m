function val=minor(A,i,j)
[d,~]=size(A);
if i==1
    if j==1
        val=A(2:d,2:d);
    else
        val=A(2:d,[1:(j-1) (j+1):d]);
    end
else
    if j==1
        val=A([1:(i-1) (i+1):d],2:d);
    else
        val=A([1:(i-1) (i+1):d],[1:(j-1) (j+1):d]);
    end
end
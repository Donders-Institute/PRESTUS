function B = changem_vectorized(A,newval,oldval)

    B = A;
    [valid,id] = max(bsxfun(@eq,A(:),oldval(:).'),[],2); %//'
    B(valid) = newval(id(valid));

return;
function La = graph_reg_update2(L,Lam,inv)
%
%
N = length(Lam);
La = cell(1,N);
for i=1:N
    La{i} = inv{i}*t2m(L-Lam{i}, i);
    s = size(L);
    modc = setdiff(1:N, i);
    La{i} = ipermute(reshape(La{i}, [s(i), s(modc)]), [i, modc]);
end

end
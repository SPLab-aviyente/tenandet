function U = uInit(sizes, modes, ranks)
%uInit Initialize U's for mapping
%  U = uInit(size, modes, ranks)
%  sizes : size of the tensor to be mapped.
%  modes : modes along which the tensor will be mapped.
%  ranks : TT-ranks of the factors.
l  = length(modes);
U  = cell(1, l);

U{1} = randn(1, sizes(modes(1)), ranks(1));
U{end} = randn(ranks(end), sizes(modes(1)), 1);
for i=2:l-1
    U{i} = randn(ranks(i-1), sizes(modes(i)), ranks(i));
end
U(l/2+2:l+1) = U(l/2+1:end);
U{l/2+1} = randn(ranks(l/2), 1, ranks(l/2));
end


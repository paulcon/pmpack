function I = index_set(type,order,dim,constraint)

if ~exist('dim','var') || isempty(dim) 
    dim=length(order);
else
    if length(order)>1 && dim~=length(order)
        error('dim is not consistent with order.');
    end
end

if isequal(type,'tensor')
    if isscalar(order), order=order*ones(dim,1); end
    I=1;
    for i=1:dim
        I=[kron(I,ones(order(i)+1,1)) kron(ones(size(I,1),1),(0:order(i))')];
    end
    I=I(:,2:end);
    I=I';
elseif isequal(type,'full') || isequal(type,'complete') || isequal(type,'total order')
    if ~isscalar(order), error('Order must be a scalar for constrained index set.'); end
    constraint=@(a) sum(a)<=order;
    I=index_set('constrained',order,dim,constraint);
elseif isequal(type,'constrained')
    if ~exist('constraint','var') || isempty(constraint), error('No constraint provided.'); end
    if ~isscalar(order), error('Order must be a scalar for constrained index set.'); end
    I=[]; 
    index=zeros(1,dim); limit=order*ones(1,dim);
    while any(index(:)<limit(:))
        if constraint(index(:))
            I=[I index(:)];
        end
        index=increment(index,limit);
    end
    if constraint(limit(:)), I=[I limit(:)]; end
else 
    error('Unrecognized type: %s',type);
end

function r = increment(index, limit)

dimension = length(index);

for i=dimension:-1:1
    if (index(i) < limit(i))
        index(i) = index(i) + 1;
        break;
    else
        index(i) = 0;
    end 
end

r = index;
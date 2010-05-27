%% Check the quadrature rules

msgid = 'pmpack:quadratureTests';

%% Test jacobi quadrature
% For Jacobi quadrature, 
for i=1:100
    % pick random jacobi paramters
    ab=floor(16*rand(2,1)); lr=sort(rand(2,1));
    s=jacobi_parameter(lr(1),lr(2),ab(1),ab(1));
    N=2*floor(16*rand(1))+1; [x,w] = gaussian_quadrature(s,N); xw = [x,w];
    midpt = xw(floor(size(xw,1)/2)+1,1);
    ex = (lr(2)-lr(1))/2+lr(1);
    if abs(midpt - ex) > 100*eps,
        error(msgid,...
            'gaussian_quadrature(%i) misses midpoint for jacobi_parameter(%i,%g,%g,%g,%g)',...
            N,lr(1),lr(2),ab(1),ab(2));
    end
    if abs(norm(xw(:,2),1)-1) > 100*eps, 
        error(msgid,...
            'gaussian_quadrature(%i,jacobi_parameter(%g,%g,%g,%g)) norm incorrect',...
            N,lr(1),lr(2),ab(1),ab(2)); 
    end
end

% function for creating animated gif of polynomial convergence in
% twobytwo_demo.m

P = twobytwo_func(0.2);
iAb = P.solve; 
s = P.s;
x1 = @(t,e) (2-t)./(1+e-t.^2);
x2 = @(t,e) (1+e-2*t)./(1+e-t.^2);
ss = linspace(-1,1,501);

X = pseudospectral(iAb,s,1);

   
x_interp = zeros(2,length(ss));
for j=1:length(ss)
    x_interp(:,j) = evaluate_expansion(X,ss(j));
end

gp = gaussian_quadrature(s,2);
subplot(1,2,2);
plot(ss,x2(ss,0.2),'b-',...
    ss,x_interp(2,:),'r-',...
    gp,x2(gp,0.2),'bo');
legend('Exact','Interp','Gauss Points','Location','NorthEast');
ylim([-5 20]);
set(gca,'nextplot','replacechildren');


subplot(1,2,1);
plot(ss,x1(ss,0.2),'b-',...
    ss,x_interp(1,:),'r-',...
    gp,x1(gp,0.2),'bo');
legend('Exact','Interp','Gauss Points','Location','NorthEast');
ylim([-5 20]);
set(gca,'nextplot','replacechildren');

set(gcf,'Color','white');

f = getframe(gcf);
[im,map]=rgb2ind(f.cdata,256,'nodither');
im(1,1,1,15)=0;



for i=1:15
    X = pseudospectral(iAb,s,i);
    
    x_interp = zeros(2,length(ss));
    for j=1:length(ss)
        x_interp(:,j) = evaluate_expansion(X,ss(j));
    end
    
    gp = gaussian_quadrature(s,i+1);
    
    subplot(1,2,2);
    plot(ss,x2(ss,0.2),'b-',...
        ss,x_interp(2,:),'r-',...
        gp,x2(gp,0.2),'bo');
    legend('Exact','Interp','Gauss Points','Location','NorthEast');
    xlabel('s'); ylabel('x_2(s)');
    ylim([-5 20]);
    
    
    subplot(1,2,1);
    plot(ss,x1(ss,0.2),'b-',...
        ss,x_interp(1,:),'r-',...
        gp,x1(gp,0.2),'bo');
    legend('Exact','Interp','Gauss Points','Location','NorthEast');
    xlabel('s'); ylabel('x_1(s)');
    ylim([-5 20]);
    
    f = getframe(gcf);
    im(:,:,1,i)=rgb2ind(f.cdata,map,'nodither');

end

imwrite(im,map,'convergence.gif','DelayTime',0.5,'LoopCount',inf);






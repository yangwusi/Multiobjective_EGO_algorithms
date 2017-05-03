function y=ZDT4(x, m)
% 0<x1<1
% -5<xi<5  i=2,3,...10
y(:,1)=x(:,1);
g=1+10*(size(x,2)-1)+sum(x(:,2:end).^2-10*cos(4*pi*x(:,2:end)),2);
h=1-(y(:,1)./g).^0.5;
y(:,2)=g.*h;
end
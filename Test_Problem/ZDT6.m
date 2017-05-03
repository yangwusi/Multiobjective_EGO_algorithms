function y=ZDT6(x, m)
% 0<xi<1

y(:,1)=1-exp(-4*x(:,1)).*(sin(6*pi*x(:,1))).^6;
g=1+9*(sum(x(:,2:end),2)/(size(x,2)-1)).^0.25;
h=1-(y(:,1)./g).^2;
y(:,2)=g.*h;
end
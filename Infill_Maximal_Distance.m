function obj = Infill_Maximal_Distance(x, sample_x)
num_point = size(x,1);
obj = zeros(num_point,1);
for ii = 1 : num_point
    obj(ii,:) = -min(sum((repmat(x(ii, :),size(sample_x),1) - sample_x).^2,2));
end

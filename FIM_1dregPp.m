function FIMpp = FIM_1dregPp(x,p0)
% fim_3dlinreg(x): return Fischer Information Matrix terms given design space x
% for a 3D linear regression model
p=p0;

Z = power(x(:), 0:p);  
   
q = size(Z,2);

FIMpp = zeros(q,q,length(x));

for i = 1:length(x)
    FIMpp(:,:,i) = Z(i,:)'*Z(i,:);
end

end
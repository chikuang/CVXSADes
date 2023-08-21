function FIMp3 = FIM_1dregP3(x)
% fim_3dlinreg(x): return Fischer Information Matrix terms given design space x
% for a 3D linear regression model
p=3;

Z = power(x(:), 0:p);  
   
q = size(Z,2);

FIMp3 = zeros(q,q,length(x));

for i = 1:length(x)
    FIMp3(:,:,i) = Z(i,:)'*Z(i,:);
end

end
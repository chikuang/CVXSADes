function FIM_polyP = FIM_polyP(x, p)
  % fim_3dlinreg(x): return Fisher Information Matrix terms given design space x
  % for a 3D linear regression model
  Z = power(x(:), 0:p);  
     
  q = size(Z,2);
  
  FIM_polyP = zeros(q,q,length(x));
  
  for i = 1:length(x)
      FIM_polyP(:,:,i) = Z(i,:)' * Z(i,:);
  end
end
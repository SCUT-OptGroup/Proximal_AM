function [Z,obj] = centroid_meancost( bdim, supp, w, c0 )

  n = length(bdim);
  
  m = length(c0.w);
  
  temp_idx = [1,cumsum(bdim)+1];
  
  D = zeros(n,1);  
  
  Z= zeros(m,sum(bdim));
  
  for i=1:n                   
      
    idx = temp_idx(i):temp_idx(i+1)-1;  
    
    [D(i), Z(:,idx)] = centroid_LP(c0.supp, c0.w, supp(:,idx), w(idx));
    
  end
  
  obj = mean(D);
  
  fprintf('\n \n Objective value:'); 
  fprintf (' %f \n',obj); 
  
end


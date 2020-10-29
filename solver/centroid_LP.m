function [fval, x] = centroid_LP(X, wX, Y, wY)

  default_options = optimset('Display','off', 'Diagnostics','off');
  
  n = size(X,2);
  
  m = size(Y,2);
  
  wX = wX/sum(wX);
  
  wY = wY/sum(wY);

  D = Fmap(X, Y);
  
  f = reshape(D, n*m, 1);

  bb = repmat({ones(1,n)},m,1);
  
  BB = blkdiag(bb{:});
  
  Aeq = [repmat(eye(n,n),1,m); BB];
  
  beq = [wX'; wY'];
  
  [x, fval] = linprog(f, [], [], Aeq, beq, zeros(n*m,1), [], [], default_options );

  x = reshape(x, n, m);
  
  x(x<0) = 0;
  
end

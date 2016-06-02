N = 20;
y = zeros(N);

zs = linspace(0,1,N);

tmp1 = @(x) smin(x,0.5,0.1);

for i=1:100
  x=randn();
  [f1,df1] = geval(tmp1,x,struct('grad_method','taylorvar'));
  [f2,df2] = tmp1(x);
  valuecheck(df1,df2);
end
  
  
for i=1:N
  for j=1:N
    zij = [zs(i);zs(j)];
    y(N-i+1,j) = smin(norm(zij-[1;0]),0.5,0.1) + smin(norm(zij-[0;1]),0.5,0.1);
  end
end

imagesc(y)
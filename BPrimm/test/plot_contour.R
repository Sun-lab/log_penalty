
x1 = seq(-0.1, 0.5, by=0.005)
y1 = seq(-0.1, 0.5, by=0.005)
x0 = 0.3
y0 = 0.2
w1 = 5
w2 = 1
lambda = 0.2

z1 = matrix(nrow=length(x1), ncol=length(y1))

for(i in 1:length(x1)){
  xi = x1[i]
  for(j in 1:length(y1)){
    yj = y1[j]
  	z1[i,j] = (xi-x0)^2/w1 + (yj-y0)^2/w2 + lambda*sqrt(xi^2 + yj^2)
  }
}

filled.contour(x=x1, y=y1, z=z1, color = terrain.colors)

points(0, 0, pch=2, col="red")
points(x0, y0, pch=3, col="blue")

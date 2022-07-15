## RANGE KUTTA METHOD
## Roll: bs2016
library("scatterplot3d");

rk = function(x0,y0,vx0,vy0){
  
n = 500;
h = 0.1;

t = rep(0,n);     
x = rep(0,n);      # x,y,z are co-ordinates of the point mass at any given time
y = rep(0,n);
z = rep(0,n);
u = rep(0,n);      # u = distance between origin and the projection of the position vector on the xy plane
vx = rep(0,n);     # x1,y1,u1 are the first derivatives wrt t
vy = rep(0,n);
vu = rep(0,n);

z1 = rep(0,n);     # first derivative of z wrt u
z2 = rep(0,n);     # second derivative of z wrt u

R = rep(0,n);      # R = normal force

g = 9.8;

x[1] = x0;         # setting initial values
y[1] = y0;

u[1] = sqrt(x[1]^2 + y[1]^2);
z[1] = sqrt(u[1]-1);
vx[1] = vx0;
vy[1] = vy0;
vu[1] = (x[1]*vx[1] + y[1]*vy[1])/u[1];
z1[1] = 1/(2*z[1]);
z2[1] = -1/(4*(z[1]^3));
R[1] = (z1[1]*(vx[1]^2 + vy[1]^2 - vu[1]^2)/u[1] + vu[1]^2*z2[1] + g)/(u[1]*(z1[1] + 1/z1[1]));

# we solve using 4th order Range Kutta method

for(i in 2:n){
        
  t[k] = t[k-1] + h;
  
  e1 = h*vx[k-1];                  
  f1 = h*vy[k-1];
  i1 = -h*R[k-1]*x[k-1];
  j1 = -h*R[k-1]*y[k-1];

  e2 = h*(vx[k-1] + i1/2);
  f2 = h*(vy[k-1] + j1/2);
  i2 = -h*R[k-1]*(x[k-1] + e1/2);
  j2 = -h*R[k-1]*(y[k-1] + f1/2);

  e3 = h*(vx[k-1] + i2/2);
  f3 = h*(vy[k-1] + j2/2);
  i3 = -h*R[k-1]*(x[k-1] + e2/2);
  j3 = -h*R[k-1]*(y[k-1] + f2/2);

  e4 = h*(vx[k-1] + i3);
  f4 = h*(vy[k-1] + j3);
  i4 = -h*R[k-1]*(x[k-1] + e3);
  j4 = -h*R[k-1]*(y[k-1] + f3);
  
  x[k] = x[k-1] + (e1+2*e2+2*e3+e4)/6;
  y[k] = y[k-1] + (f1+2*f2+2*f3+f4)/6;
  vx[k] = vx[k-1] + (i1+2*i2+2*i3+i4)/6;
  vy[k] = vy[k-1] + (j1+2*j2+2*j3+j4)/6;
  u[k] = sqrt(x[k]*x[k] + y[k]*y[k]);
  z[k] = sqrt(u[k]-1);
  vu[k] = (x[k]*vx[k] + y[k]*vy[k])/u[k];
  
  z1[k] = 1/(2*z[k]);
  z2[k] = -1/(4*z[k]*z[k]*z[k]);
  
  R[k] = (z1[k]*(vx[k]^2 + vy[k]^2 - vu[k]^2)/u[k] + vu[k]^2*z2[k] + g)/(u[k]*(z1[k] + 1/z1[k]));
}

# creating the grid

a = seq(-3.25,3.25,len = 1000);

a1 = sqrt(10.5625 - a^2);
a2 = rep(1.5,1000);

b = seq(-18,18,len = 1000);

b1 = sqrt(324 - b^2);
b2 = rep(sqrt(17),1000);

c = seq(1.5,sqrt(17),len = 1000);

c1 = 1 + c^2;
c2 = c1/sqrt(2);
c3 = rep(0,1000);

# plotting the graph in 3D

for(k in 1:n){
s <- scatterplot3d(a,a1,a2,grid = FALSE,axis = FALSE,angle = 85,
                   zlim = c(0,5),xlim = c(-18,18),ylim = c(-18,18),ty = 'l')
s$points3d(a,-a1,a2,ty = 'l')
s$points3d(b,b1,b2,ty = 'l')       
s$points3d(b,-b1,b2,ty = 'l')
s$points3d(c1,c3,c,ty = 'l')
s$points3d(-c1,c3,c,ty = 'l')
s$points3d(c3,c1,c,ty = 'l')
s$points3d(c3,-c1,c,ty = 'l')
s$points3d(c2,c2,c,ty = 'l')       
s$points3d(-c2,-c2,c,ty = 'l')
s$points3d(c2,-c2,c,ty = 'l')
s$points3d(-c2,c2,c,ty = 'l')

# plot the point mass
s$points3d(x[k],y[k],z[k],pch = 16,cex = 2,col = "blue")
Sys.sleep(0.000001)
}
# list of time and position vectors
return(list(time = t,x = x,y = y,z = z))    
}

# for the initial conditions given
result = runge_kutta(10,0,0,5)



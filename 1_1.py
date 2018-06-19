#encoding=utf-8
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
from math import pi
from math import pow
import cmath
from math import cos
from math import sin
from math import e
from math import atan



mesh=Mesh("mesh.xml")
def u0_boundary(x):
	return near(x[1], 50000) 
# no multiply i
def sesq_inner(u_r, u_i, v_r, v_i):

	re_part = inner(u_r, v_r) - inner(u_i, v_i)
	im_part =  inner(u_r, v_i) + inner(u_i, v_r)	
	return	re_part + im_part
#multiply i
def sesq_inner1(u_r, u_i, v_r, v_i):
    re_part = -inner(u_r, v_i) - inner(u_i, v_r)
    im_part =  inner(u_r, v_r) - inner(u_i, v_i)
    return re_part + im_part

class Bottom(SubDomain):
		def inside(self,x,on_boundary):
			return near(x[1],-100000.0,DOLFIN_EPS) 
bottom=Bottom()

class Omega0(SubDomain):
	def inside(self, x, on_boundary):
		return x[1]>=0.0 
class Omega1(SubDomain):
	def inside(self, x, on_boundary):
		return 2*x[1]+x[0]>=0.0 and x[1]>=-4000 and x[1]<=-1000
class Omega2(SubDomain):
	def inside(self, x, on_boundary):
		return x[1]>=-1000.0 and x[1]<=0.0 and x[0]>=-100000.0 and x[0]<=100000.0 

boundaries=FacetFunction("size_t",mesh,1)
bottom.mark(boundaries,2)

sigma0=1.0/pow(10,10)
sigma1=1.0/30
sigma2=1.0/600


e1=8.85*pow(10,-12)
miu=4.0*pi*pow(10,-7)
p_a=[]
w1=[]
u0=Constant((1.0,0))
f=open("aaa.txt","w")
f1=open("bbb.txt","w")

ds = Measure('ds', domain=mesh, subdomain_data=boundaries)
dx = Measure('dx', domain=mesh, subdomain_data=subdomains)



for i in range (-20,21,1):
	i=i*0.1
	w=2*pi/pow(10,i)
	w1.append(i)
	print(w)
	V= FunctionSpace(mesh, "CG",	2)
	W=V*V



	u_r, u_i =TrialFunction(W)
	v_r, v_i =TestFunction(W)
	
	a= (1/(w*miu))*sesq_inner1(grad(u_r), grad(u_i), grad(v_r), grad(v_i))*dx-miu*e1*sesq_inner1(u_r, u_i, v_r, v_i)*dx
	
	l =sigma0*sesq_inner(u_r, u_i, v_r, v_i)*dx(0)+sigma1*sesq_inner(u_r, u_i, v_r, v_i)*dx(1)+sigma2*sesq_inner(u_r, u_i, v_r, v_i)*dx(2)
	c=(sqrt(sigma2)/sqrt(2*w*miu))*(sesq_inner1(u_r, u_i, v_r, v_i)-sesq_inner(u_r, u_i, v_r, v_i))*ds(2)
 



	bcs0 = DirichletBC(W,u0,u0_boundary)  # Real BC


	F=a-l+c
	A=lhs(F)
	L=rhs(F) 
	u = Function(W)
	solve(A == L, u, bcs0)
	u_r,u_i=u.split(u)

	
	W= VectorFunctionSpace(mesh, 'Lagrange', 1)
	grad_u_r= project(grad(u_r), W)
	grad_u_i= project(grad(u_i), W)
	for j in range(0,12100,100):
		grad_u=complex(grad_u_r(j,0)[1],grad_u_i(j,0)[1])
		u=complex(u_r(j,0),u_i(j,0))
		p=w*miu*abs(u/grad_u*u/grad_u)
		p_a.append(p)
		Z_TE=u*complex(0,w*miu)/grad_u
		sita_TE=180/pi*atan(Z_TE.imag/Z_TE.real)
		f.write(str(-i)+"\t"+str(j/1000.0)+"\t"+str(p)+"\n")
		f1.write(str(-i)+"\t"+str(j/1000.0)+"\t"+str(sita_TE)+"\n")
f.close
f1.close
file=File("subdomain.pvd")
file<<subdomains

interactive()

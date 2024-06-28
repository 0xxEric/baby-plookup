import random
import galois
import numpy as np
from fft import fft
from Eric_poly_commit import PolyCommitment
from poly_utils import PrimeField
from typing import Any

#region Step1: Initial parameter settings 初始参数设置
# α β γ ζ等随机参数：均于交互过程中设置，不在初始参数中设置。
modulus = 101 #取一个较小的域，来进行计算 
pf=PrimeField(modulus)
pc=PolyCommitment(modulus)
GFM=galois.GF(modulus) #指的是以modulus的类
g=7
G = pc.pf.exp(g, (pc.pf.modulus-1) // 4)

def lagrange_poly_in_finite_field(x, y):
    return galois.lagrange_poly(GFM(x), GFM(y))
def get_div_num(a, b):
    return GFM(a) / GFM(b)

class plonkup:
    def __init__(self,modulus):
        self.modulus = modulus
        self.pf=PrimeField(modulus)
        self.pc=PolyCommitment(modulus)
        self.g=7
        self.G = pc.pf.exp(g, (pc.pf.modulus-1) // 4)
        self.beta=GFM(random.randint(0, modulus - 1))
        self.gama=GFM(random.randint(0, modulus - 1))

    def Gab(self,a,b):
        return a+self.beta*b+self.gama*(GFM(1)+self.beta)
    
    def shift(self,t):
        l=len(t)
        m=[]
        for i in range(l-1):
          m.append(t[i+1])
        m.append(t[0])
        return m
    
pl=plonkup(modulus)
print(f"beta2={pl.beta.max()},gama2={pl.gama.max()}")
roots_of_unity = GFM([1, 10, 100, 91]) #真实项目中，应该根据约束、向量t及f的长度，来确定degree、其roots的个数，本处仅为demo示例，直接给出roots
#endregion


#region step2 input vectors: t,f,s,s_even,s_odd
#这里就直接使用示例，把各个向量写出来了。真实项目中应该生成出来
t=GFM([1,2,3,4])
f=GFM([2,4,4,1])
s=GFM([1,1,2,2,3,4,4,4])
s_even=GFM([1,2,3,4])
s_odd=GFM([1,2,4,4])
F=GFM([[2,2],[4,4],[4,4],[1,1]])
T=GFM([[1,2],[2,3],[3,4],[4,1]])
Seven=GFM([[1,1],[2,2],[3,4],[4,4]])
Sodd=GFM([[1,2],[2,3],[4,4],[4,1]])
#endregion


#region step3 calculate the values and coeffs
#Calculate  values f_values and g_values 
f_values = []
g_values = []
q_values = []
for i in range(len(roots_of_unity)):
    f_value=pl.Gab(F[i][0],F[i][1])*pl.Gab(T[i][0],T[i][1])    
    g_value=pl.Gab(Seven[i][0],Seven[i][1])*pl.Gab(Sodd[i][0],Sodd[i][1])
    f_values.append(f_value)
    g_values.append(g_value)
    q_values.append(get_div_num(f_value,g_value))

# Calculate accumulator values z_values
z_values = [1]
for i in range(len(roots_of_unity) - 1):
    z_values.append(z_values[-1] * get_div_num(f_values[i], g_values[i]))
assert z_values[len(roots_of_unity)-1]*q_values[len(roots_of_unity)-1]==1 #验证zN=1

# Calculate z_w_values by shift z_values, the last value is 1
z_w_values = z_values[1:] + z_values[:1]

# other values
L_0_values = GFM([1, 0, 0, 0])
t_w_values=pl.shift(t) #将t进行移位，形成t_w
seven_w_values=pl.shift(s_even) #将seven进行移位，形成s_even
# construct the polys
f_poly = lagrange_poly_in_finite_field(roots_of_unity, f) 
t_poly = lagrange_poly_in_finite_field(roots_of_unity, t) 
t_w_poly=lagrange_poly_in_finite_field(roots_of_unity, t_w_values) 
s_even_poly = lagrange_poly_in_finite_field(roots_of_unity, s_even) 
s_w_even_poly = lagrange_poly_in_finite_field(roots_of_unity, seven_w_values) 
s_odd_poly = lagrange_poly_in_finite_field(roots_of_unity, s_odd) 
z_poly=lagrange_poly_in_finite_field(roots_of_unity, z_values) 
z_w_poly=lagrange_poly_in_finite_field(roots_of_unity, z_w_values) 
L_0_poly = lagrange_poly_in_finite_field(roots_of_unity, L_0_values)
g3_poly=(GFM(1)+pl.beta)*(f_poly+pl.gama)*(t_poly+pl.beta*t_w_poly+pl.gama*(GFM(1)+pl.beta))
g4_poly=(s_even_poly+pl.beta*s_odd_poly+pl.gama*(GFM(1)+pl.beta))*(s_odd_poly+pl.beta*s_w_even_poly+pl.gama*(GFM(1)+pl.beta))
#endregion


#region step4 calculate C1\C2\combined poly，and check
#Calculate  values f_values and g_values 
C1_poly=L_0_poly * (z_poly - GFM(1))  #即L0(x)*(Z(x)-1)
C2_poly=z_w_poly*g4_poly-z_poly*g3_poly 

alpha2=GFM(random.randint(0, modulus - 1))
print(f"alpha2={alpha2.max()}")
combined_poly=alpha2*C1_poly+C2_poly #使用alpha聚合两式
print("Vanishing polynomial combined_poly: ", combined_poly)

# Vanishing polynomial: (X - 1)(X - 10)(X - 100)(X - 91) #也可以直接用X^N-1 ?
z_H_poly = galois.Poly([1, roots_of_unity[0]], field = GFM) * galois.Poly([1, roots_of_unity[1]], field = GFM) * galois.Poly([1, roots_of_unity[2]], field = GFM) * galois.Poly([1, roots_of_unity[3]], field = GFM)
print("Vanishing polynomial z_H_poly: ", z_H_poly)

quot_poly = combined_poly // z_H_poly #计算商多项式
print("Quotient polynomial quot_poly: ", quot_poly)
assert combined_poly == quot_poly * z_H_poly, "Division has a remainder(should not have)"
print("successful checked!")
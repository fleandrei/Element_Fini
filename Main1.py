# -*- coding:Utf-8 -*
import gmsh
import sys
import numpy as np
from Maillage import *
from Assemblage_Matrice import *
from Functions import *
from scipy.sparse import linalg
import matplotlib.pyplot as plt


gmsh.initialize(sys.argv)
h=0.25
gener_Gmsh_carre("Carre_Unite.msh", h)#On génére le maillage Carré_Unité.msh


Tag_Espace=2 #Tag physique de l'espace sur lequel est définit la fonction
Tag_Dirichlet=[0,1] #Tag physique du domaine de dirichlet

Mat=Triplet() #Matrice A sous forme de triplet
Msh=Mesh()
Msh.GmshToMesh("Carre_Unite.msh")

nbr_Points=len(Msh.Points)
#Ones=[1 for i in range(nbr_Points)]

Mass(Msh, 2, Tag_Espace, Mat)
#print(T)

Mat_Mass=(sparse.coo_matrix((Mat.data[0],(Mat.data[1][0],Mat.data[1][1])))).tocsr()
Mat_Mass=Mat_Mass[1:nbr_Points+1,1:nbr_Points+1]
#print(sum(Mat_Mass*Ones))


Rigid(Msh, Tag_Espace, Mat)
#print(T)

Mat_Rigi=(sparse.coo_matrix(Mat.data)).tocsr()
Mat_Rigi=Mat_Rigi[1:nbr_Points+1,1:nbr_Points+1]
#print(Mat_Rigi*Ones)
T2=Triplet()
Rigid(Msh, Tag_Espace, T2)
Mat_Rigi=(sparse.coo_matrix(T2.data)).tocsr()
Mat_Rigi=Mat_Rigi[1:nbr_Points+1,1:nbr_Points+1]

#print((30 in T.data[1][0]) or (30 in T.data[1][1]))
B=Gener_B(Msh,Tag_Espace,F,2)
Dirichlet(Msh,Tag_Dirichlet,G, Mat,B)

A=(sparse.coo_matrix(Mat.data)).tocsr()
A=A[1:nbr_Points+1,1:nbr_Points+1]
#print(A.toarray())
U=sparse.linalg.spsolve(A, B) #Approximation de la solution

Ainv=sparse.linalg.spsolve(A, np.eye(nbr_Points))

indtest=18
print("Matrice mass\n"+str(Mat_Mass.toarray()[indtest]))
print("Matrice rigi\n"+str(Mat_Rigi.toarray()[indtest]))
print("Matrice A\n"+str(A.toarray()[indtest]))
print("Matrice A inv\n"+str(Ainv[indtest]))
print(U)
print(B)
print(A*U)
#print(B)

X=[p.X for p in Msh.Points]
Y=[p.Y for p in Msh.Points]
#print(X)
#print(Y)


#print(F_discret)


connectivity=[]
for tri in Msh.Triangles:
  connectivity.append([ p.ID_Glob -1 for p in tri.Points]) 
#print(connectivity)

plt.tricontourf(X, Y, connectivity, U, 17)
plt.colorbar()
plt.show()


Uref = np.zeros((nbr_Points,))
for pt in Msh.Points:
  I = int(pt.ID_Glob)
  Uref[I -1] = Solution(pt.X, pt.Y)

print(Uref)
Bref=A*Uref
print(Uref-U)

plt.figure()
plt.tricontourf(X, Y, connectivity, Uref, 17)
plt.colorbar()
plt.show()
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
gener_Gmsh_carre(h, "Carre_Unite.msh")#On génére le maillage Carré_Unité.msh


Tag_Espace=2 #Tag physique de l'espace sur lequel est définit la fonction
Tag_Dirichlet=[0,1] #Tag physique du domaine de dirichlet

Mat=Triplet() #Matrice A sous forme de triplet
Msh=Mesh()
Msh.GmshToMesh("Carre_Unite.msh")

nbr_Points=len(Msh.Points)

Mass(Msh, 2, Tag_Espace, Mat)

Mat_Mass=(sparse.coo_matrix((Mat.data[0],(Mat.data[1][0],Mat.data[1][1])))).tocsr()
Mat_Mass=Mat_Mass[1:nbr_Points+1,1:nbr_Points+1]


Rigid(Msh, Tag_Espace, Mat)

Mat_Rigi=(sparse.coo_matrix(Mat.data)).tocsr()
Mat_Rigi=Mat_Rigi[1:nbr_Points+1,1:nbr_Points+1]

T2=Triplet()
Rigid(Msh, Tag_Espace, T2)
Mat_Rigi=(sparse.coo_matrix(T2.data)).tocsr()
Mat_Rigi=Mat_Rigi[1:nbr_Points+1,1:nbr_Points+1]

B=Gener_B(Msh,Tag_Espace,F,2)
Dirichlet(Msh,Tag_Dirichlet,G, Mat,B)

A=(sparse.coo_matrix(Mat.data)).tocsr()
A=A[1:nbr_Points+1,1:nbr_Points+1]
U=sparse.linalg.spsolve(A, B) #Approximation de la solution

Ainv=sparse.linalg.spsolve(A, np.eye(nbr_Points))

indtest=18


X=[p.X for p in Msh.Points]
Y=[p.Y for p in Msh.Points]





connectivity=[]
for tri in Msh.Triangles:
  connectivity.append([ p.ID_Glob -1 for p in tri.Points]) 

plt.tricontourf(X, Y, connectivity, U, 17)
plt.title("Solution approchee de U")
plt.colorbar()
plt.savefig("Solution_Approchée", format="png")
plt.show()


Uref = np.zeros((nbr_Points,))
for pt in Msh.Points:
  I = int(pt.ID_Glob)
  Uref[I -1] = Solution(pt.X, pt.Y)

Bref=A*Uref


plt.figure()
plt.title("Interpolation de la solution exacte")
plt.tricontourf(X, Y, connectivity, Uref, 17)
plt.colorbar()
plt.show()
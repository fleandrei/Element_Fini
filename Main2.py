# -*- coding:Utf-8 -*
import gmsh
import sys
import numpy as np
from Maillage import *
from Assemblage_Matrice import *
from Functions import *
from scipy.sparse import linalg
import matplotlib.pyplot as plt
import json


def return_Erreur(h): #Fonction permettant de renvoyer l'erreur en norme L2 
	gmsh.initialize(sys.argv)

	gener_Gmsh_carre(h,"Carre_Unite.msh")#On génére le maillage Carré_Unité.msh


	Tag_Espace=2 #Tag physique de l'espace sur lequel est définit la fonction
	Tag_Dirichlet=[0,1] #Tag physique du domaine de dirichlet

	Mat=Triplet() #Matrice A sous forme de triplet
	Msh=Mesh()
	Msh.GmshToMesh("Carre_Unite.msh")
	nbr_Points=len(Msh.Points)

	Mass(Msh, 2, Tag_Espace, Mat)
	Rigid(Msh, Tag_Espace, Mat)


	B=Gener_B(Msh,Tag_Espace,F,2)
	Dirichlet(Msh,Tag_Dirichlet,G, Mat,B)

	A=(sparse.coo_matrix(Mat.data)).tocsr()
	A=A[1:nbr_Points+1,1:nbr_Points+1]
	#print(A.toarray())
	U=sparse.linalg.spsolve(A, B) #Approximation de la solution


	X=[p.X for p in Msh.Points]
	Y=[p.Y for p in Msh.Points]


	Uref = np.zeros((nbr_Points,))
	for pt in Msh.Points:
	  I = int(pt.ID_Glob)
	  Uref[I -1] = Solution(pt.X, pt.Y)

	Point.identifiant=0
	Segment.identifiant=0
	Triangle.identifiant=0
	
	Erreur=np.linalg.norm(Uref-U,2)
	return Erreur





if len(sys.argv)==4: #Si l'utilisateur a donné 3 arguments: Hmin Hmax pas
	Hini=float(sys.argv[1])
	Hfin=float(sys.argv[2])
	Hpas=float(sys.argv[3])


	H_courant=Hini
	Erreur_log=[]
	H_log=[]
	Erreur=[]
	H=[]
	while H_courant<=Hfin:
		err=return_Erreur(H_courant)
		Erreur.append(err)
		H.append(H_courant)

		Erreur_log.append(math.log(err))
		H_log.append(math.log(H_courant))

		H_courant+=Hpas

	plt.plot(H_log, Erreur_log)
	plt.title("Evolution de l'erreur en norme L2 en fonction de h")
	plt.xlabel("ln(h)")
	plt.ylabel("ln(Erreur)")
	plt.savefig("Evolution_Erreur", format="png")
	plt.show()

	Fichier_resultat=open("Résultat_Erreur","w")
	json_string=json.dumps({ "Log de h":H_log})
	json_string=json_string[0:-1]+",\n \n"+json.dumps({"Log de l'erreur" : Erreur_log})[1:]
	json_string=json_string[0:-1]+",\n\n"+json.dumps({"h" : H})[1:]
	json_string=json_string[0:-1]+",\n\n"+json.dumps({"Erreur":Erreur})[1:]

	Fichier_resultat.write(json_string)

else:
	if len(sys.argv)!=2:
		print("Entrez la valeur du pas:\n -> "+str(sys.argv[0])+" h\n -> "+str(sys.argv[0])+" h_initial h_final pas\nValeur par défaut: h=0.25")
		h=0.25
	else:
		h=float(sys.argv[1])


	gmsh.initialize(sys.argv)

	gener_Gmsh_carre(h,"Carre_Unite.msh")#On génére le maillage Carré_Unité.msh


	Tag_Espace=2 #Tag physique de l'espace sur lequel est définit la fonction
	Tag_Dirichlet=[0,1] #Tag physique du domaine de dirichlet

	Mat=Triplet() #Matrice A sous forme de triplet
	Msh=Mesh()
	Msh.GmshToMesh("Carre_Unite.msh")

	nbr_Points=len(Msh.Points)
	Ones=[1 for i in range(nbr_Points)]

	Mass(Msh, 2, Tag_Espace, Mat)

#######Pour vérifier la matrice de masse
	Mat_Mass=(sparse.coo_matrix((Mat.data[0],(Mat.data[1][0],Mat.data[1][1])))).tocsr()
	Mat_Mass=Mat_Mass[1:nbr_Points+1,1:nbr_Points+1]
	#print(sum(Mat_Mass*Ones))


	Rigid(Msh, Tag_Espace, Mat)
	#print(T)

	Mat_Rigi=(sparse.coo_matrix(Mat.data)).tocsr()
	Mat_Rigi=Mat_Rigi[1:nbr_Points+1,1:nbr_Points+1]
	
#######Pour vérifier la matrice de Rigidité	
	T2=Triplet()
	Rigid(Msh, Tag_Espace, T2)
	Mat_Rigi=(sparse.coo_matrix(T2.data)).tocsr()
	Mat_Rigi=Mat_Rigi[1:nbr_Points+1,1:nbr_Points+1]
	#print(max(abs(Mat_Rigi*Ones)))
	#print((30 in T.data[1][0]) or (30 in T.data[1][1]))
	

	B=Gener_B(Msh,Tag_Espace,F,2)#Calcul de B
	Dirichlet(Msh,Tag_Dirichlet,G, Mat,B) #Condition de Dirichlet

	A=(sparse.coo_matrix(Mat.data)).tocsr()
	A=A[1:nbr_Points+1,1:nbr_Points+1]

	U=sparse.linalg.spsolve(A, B) #Approximation de la solution

	Ainv=sparse.linalg.spsolve(A, np.eye(nbr_Points))

	indtest=40
	#print("Matrice mass\n"+str(Mat_Mass.toarray()[indtest]))
	#print("Matrice rigi\n"+str(Mat_Rigi.toarray()[indtest]))
	#print("Matrice A\n"+str(A.toarray()[indtest]))
	#print("Matrice A inv\n"+str(Ainv[indtest]/20.7))
	#print(U)
	#print(B)

	X=[p.X for p in Msh.Points]
	Y=[p.Y for p in Msh.Points]
	

	connectivity=[]
	for tri in Msh.Triangles:
	  connectivity.append([ p.ID_Glob -1 for p in tri.Points]) 

	fig_sol_approx=plt.tricontourf(X, Y, connectivity, U, 12)

	plt.title("Solution approchee de U")
	plt.colorbar()
	plt.show()


	Uref = np.zeros((nbr_Points,))
	for pt in Msh.Points:
	  I = int(pt.ID_Glob)
	  Uref[I -1] = Solution(pt.X, pt.Y)

	#print(Uref)
	Bref=A*Uref
	#print(Uref-U)


	Erreur=np.linalg.norm(Uref-U,2)

	print("L'erreur L2 avec un h de {} est de {}".format(h,Erreur))
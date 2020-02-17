# -*- coding:Utf-8 -*
import gmsh
import sys
import numpy as np

from Maillage import *

def mass_elem(element, triplets, alpha=1.): #Calcule la matrice de masse élémentaire
	jac=element.jac()
	JacsurDouze=jac/12
	JacsurVingt4=jac/24
	
	if element.Name=="Segment":
		I=element.Points[0].ID_Glob
		J=element.Points[1].ID_Glob
		triplets.append(I,I, JacsurDouze)
		triplets.append(I,J, JacsurVingt4)

		triplets.append(J,I, JacsurVingt4)
		triplets.append(J,J, JacsurDouze)

	if element.Name=="Triangle":
		I=element.Points[0].ID_Glob
		J=element.Points[1].ID_Glob
		K=element.Points[2].ID_Glob

		triplets.append(I,I, JacsurDouze)
		triplets.append(I,J, JacsurVingt4)
		triplets.append(I,K, JacsurVingt4)

		triplets.append(J,I, JacsurVingt4)
		triplets.append(J,J, JacsurDouze)
		triplets.append(J,K, JacsurVingt4)

		triplets.append(K,I, JacsurVingt4)
		triplets.append(K,J, JacsurVingt4)
		triplets.append(K,K, JacsurDouze)

	return triplets


def Mass(msh, dim, physical_tag, triplets):
	Elements=msh.getElements(dim, physical_tag)
	nbrElem=len(Elements)
	for i in range(nbrElem):
		triplets=mass_elem(Elements[i], triplets)

	return triplets


def Jacobien(Elem):
	jac=np.array([[Elem.Points[1].X - Elem.Points[0].X, Elem.Points[2].X - Elem.Points[0].X],[Elem.Points[1].Y - Elem.Points[0].Y, Elem.Points[2].Y - Elem.Points[0].Y]])
	return jac


def Bp(Elem, jacob=None):#Calcule la matrice de passage Bp 
	if jacob==None:
		jacob=Jacobien(Elem)

	detjac=jacob[0][0]*jacob[1][1] - jacob[1][0]*jacob[0][1]


	BP=np.array([[Elem.Points[2].Y - Elem.Points[0].Y, Elem.Points[0].Y - Elem.Points[1].Y], [Elem.Points[0].X - Elem.Points[2].X, Elem.Points[1].X - Elem.Points[0].X]])
	BP/detjac
	return BP





def gradPhi_Chap(i): #Calcule le gradient des fonctions de forme du repère parmètrique
	if i==1:
		return np.array([-1, -1])
	elif i==2:
		return np.array([1,0])
	else:
		return np.array([0,1])

def gradPhi(element, i): #Calcule le gradient des fonctions de forme du repère normale
	BP=Bp(element)

	gradPhiChap=gradPhi_Chap(i)
	return BP.dot(gradPhiChap)



def rigi_elem(element, triplets): #Calcule la matrice rigide élémentaire
	area=element.area()
	BP=Bp(element)
	BP_prim=BP.transpose()
	Phi_chap=[gradPhi_Chap(1), gradPhi_Chap(2), gradPhi_Chap(3)]

	I=element.Points[0].ID_Glob
	J=element.Points[1].ID_Glob
	Z=element.Points[2].ID_Glob

	K=np.zeros((3,3))
	for i in range(1,4):
		for j in range(1,4):
			K[i-1,j-1]=np.dot(np.dot(Phi_chap[i-1].transpose(), np.dot(BP_prim, BP)), Phi_chap[j-1])
	K=K*area
	
	triplets.append(I,I, K[0,0])
	triplets.append(I,J, K[0,1])
	triplets.append(I,Z, K[0,2])

	triplets.append(J,I, K[1,0])
	triplets.append(J,J, K[1,1])
	triplets.append(J,Z, K[1,2])

	triplets.append(Z,I, K[2,0])
	triplets.append(Z,J, K[2,1])
	triplets.append(Z,Z, K[2,2])
	return triplets


def Rigid(msh, physical_tag, triplets):
	Elements=msh.getElements(2,physical_tag)
	nbrElem=len(Elements)
	#print(nbrElem)
	for i in range(nbrElem):
		triplets=rigi_elem(Elements[i], triplets)

	return triplets


def Gener_B(msh, physical_tag, f, order=2): #Génènre le vecteur B
	Elements=msh.getElements(2,physical_tag)
	Points=msh.getPoints(2,physical_tag)
	nbr_elem=len(Elements)
	nbr_points=len(Points)
	B=np.zeros(nbr_points)

	for i in range(nbr_elem):#On parcour les triangles
		Quad_Point=Elements[i].gaussPoint(order) #[[Poids],[(coord param 1, coord param 2)],[[coord geom x, coord geom y, coord geom z]] ]
		jac=Elements[i].jac()
		for j in range(3): #On parcourt les 3 fonctions de fomres du triangle
			contrib=0 #Contribution élémentaire
			Phiref=phiRef(j+1, Quad_Point[1]) #Valeurs prises par la fonction de forme j pour les différents points de quadratures
			for k in range(len(Quad_Point[0])):# On somme sur les points de quadratures
				contrib+=Quad_Point[0][k]*f(Quad_Point[2][k][0], Quad_Point[2][k][1])*Phiref[k]

			I=Elements[i].Points[j].ID_Glob # I=Loc2Glob(i,j)
			B[I-1]+=contrib*jac

	return B





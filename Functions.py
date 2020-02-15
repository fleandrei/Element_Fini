# -*- coding:Utf-8 -*
from Maillage import *
import numpy as np




def F(x,y):
	return (1+2*math.pi*math.pi)*math.sin(x*math.pi)*math.sin(y*math.pi)

def G(x,y):
	return 0

def Solution(x,y):
	return math.sin(x*math.pi)*math.sin(y*math.pi)

def Dirichlet(msh, physical_tag, g, triplets, B):
	Dirichlet_Domaine=[]
	tags=[]
	for t in physical_tag: # Le dommaine de dirichlet peut être constitué de plusieurs physical
		diri=msh.getPoints(1,t)
		lendiri=len(diri)
		i=0
		while i<lendiri:#On fait en sorte que les points qui appartiennent à plusiseurs physical ne soient pas mis en doublons
			if diri[i].ID_Glob in tags:
				diri.pop(i)
				lendiri-=1
			else:
				tags.append(diri[i].ID_Glob)
				i+=1
		Dirichlet_Domaine=Dirichlet_Domaine + diri#On récupère les points sur le bord. Ici physical_tag fait référence au dommaine de Dirichlet
	

	for i in range(len(Dirichlet_Domaine)):
		tag=Dirichlet_Domaine[i].ID_Glob
		print(tag)
		nbr_contrib=len(triplets.data[0])
		for j in range(nbr_contrib):
			if tag==triplets.data[1][0][j]:
				triplets.data[0][j]=0

		triplets.append(tag,tag,1)
		B[tag-1]=G(Dirichlet_Domaine[i].X, Dirichlet_Domaine[i].Y)
		
	nbr_contrib=len(triplets.data[0])
	j=0
	while j <nbr_contrib:#On supprime les coef qui sont égal à 0
		if triplets.data[0][j]==0:
			triplets.data[0].pop(j)
			triplets.data[1][0].pop(j)
			triplets.data[1][1].pop(j)
			nbr_contrib-=1
		else:
			j+=1

	return triplets
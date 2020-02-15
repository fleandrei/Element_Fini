# -*- coding:Utf-8 -*

import copy
import math
import gmsh
import sys
import numpy as np
from scipy import sparse

model = gmsh.model
factory = model.occ

gmsh.initialize(sys.argv)


def gener_Gmsh_carre(msh_name,h, show=False):
	model = gmsh.model
	factory = model.occ
	model.add(msh_name)
	points = []
	points.append(model.geo.addPoint(0,0,0, h,0)) #x , y, z, h, tag
	points.append(model.geo.addPoint(1,0,0, h,1))
	points.append(model.geo.addPoint(1,1,0, h,2))
	points.append(model.geo.addPoint(0,1,0, h,3))

	lines = []
	lines.append(model.geo.addLine(0, 1, 0))# point1, point 2, tag
	lines.append(model.geo.addLine(1, 2, 1))
	lines.append(model.geo.addLine(2, 3, 2))
	lines.append(model.geo.addLine(3, 0, 3))

	curveloop = model.geo.addCurveLoop([0,1,2,3])
	disk = model.geo.addPlaneSurface([curveloop])

	gmsh.model.addPhysicalGroup(1, lines, 1)
	gmsh.model.addPhysicalGroup(2, [disk], 2)

	gmsh.model.geo.synchronize()

	model.mesh.generate(2)
	gmsh.write(msh_name)

	if show:
		gmsh.fltk.run()
	
	gmsh.finalize()

	

class Triplet:
  def __init__(self):
    self.data = ([], ([], []))

  def __str__(self):
    return str(self.data)

  def append(self, I, J, val):
  	self.data[0].append(val)
  	self.data[1][0].append(I)
  	self.data[1][1].append(J)



class Point:
	identifiant=0

	def __init__(self, ID_Glob, X,Y,Z, iD=True):
		if iD:
			self.id=Point.identifiant
			Point.identifiant+=1
		else:
			self.id=None

		self.ID_Glob=ID_Glob
		self.X=X
		self.Y=Y
		self.Z=Z

	def __str__(self):
		return "Global ID= "+str(self.ID_Glob)+", id= "+str(self.id)+", X="+str(self.X)+", Y="+str(self.Y)+", Z="+str(self.Z)+"\n"


class Segment:
	identifiant=0

	def __init__(self, Points, Tag_Phys, id=True):
		if id:
			self.id=Segment.identifiant
			Segment.identifiant+=1
		else:
			self.id=None

		self.Tag_Phys=Tag_Phys
		self.Points=[]
		self.Points.append(copy.copy(Points[0]))
		self.Points.append(copy.copy(Points[1]))
		self.Aire=None
		self.name="Segment"


	def area(self):
		if self.Aire==None:
			self.Aire=math.sqrt(math.pow((self.Points[0].X - self.Points[1].X),2) + math.pow((self.Points[0].Y - self.Points[1].Y), 2) + math.pow((self.Points[0].Z - self.Points[1].Z), 2) )
		
		return self.Aire

	def jac(self):
		return self.area()

class Triangle:
	identifiant=0

	def __init__(self, Points, Tag_Phys, id=True):
		if id:
			self.id=Triangle.identifiant
			Triangle.identifiant+=1
		else:
			self.id=None
		
		self.Tag_Phys=Tag_Phys
		self.Points=[]
		self.Points.append(copy.copy(Points[0]))
		self.Points.append(copy.copy(Points[1]))
		self.Points.append(copy.copy(Points[2]))
		self.Aire=None
		self.Name="Triangle"

	def area(self):
		if self.Aire==None:
			S1=Segment(self.Points[0:2], 0, False)
			S2=Segment(self.Points[1:], 0, False)
			S3=Segment([self.Points[0], self.Points[2]], 0, False)
			a=S1.area()
			b=S2.area()
			c=S3.area()
			s=0.5*(a+b+c)
			self.Aire=math.sqrt(s*(s-a)*(s-b)*(s-c))

		return self.Aire

	def jac(self):
		return 2*self.area()

	

	def gaussPoint(self, order=2):
		Param=[[],[],[]]
		
		if order ==1:
			Param[0].append(1./6)
			Param[1].append((1./3, 1./3))
			Phi1=Phi_chap(1, 1./3, 1./3)
			Phi2=Phi_chap(2, 1./3, 1./3)
			Phi3=Phi_chap(3, 1./3, 1./3)
			X=Phi1*self.Points[0].X + Phi2*self.Points[1].X + Phi3*self.Points[2].X
			Y=Phi1*self.Points[0].Y + Phi2*self.Points[1].Y + Phi3*self.Points[2].Y
			Z=Phi1*self.Points[0].Z + Phi2*self.Points[1].Z + Phi3*self.Points[2].Z
			Param[2].append([X,Y,Z])

		if order==2:

			Param[0].append(1./6)
			Param[1].append((1./6, 1./6))
			Phi1=Phi_chap(1, 1./6, 1./6)
			Phi2=Phi_chap(2, 1./6, 1./6)
			Phi3=Phi_chap(3, 1./6, 1./6)
			X=Phi1*self.Points[0].X + Phi2*self.Points[1].X + Phi3*self.Points[2].X
			Y=Phi1*self.Points[0].Y + Phi2*self.Points[1].Y + Phi3*self.Points[2].Y
			Z=Phi1*self.Points[0].Z + Phi2*self.Points[1].Z + Phi3*self.Points[2].Z
			Param[2].append([X,Y,Z])

			Param[0].append(1./6)
			Param[1].append((4./6, 1./6))
			Phi1=Phi_chap(1, 4./6, 1./6)
			Phi2=Phi_chap(2, 4./6, 1./6)
			Phi3=Phi_chap(3, 4./6, 1./6)
			X=Phi1*self.Points[0].X + Phi2*self.Points[1].X + Phi3*self.Points[2].X
			Y=Phi1*self.Points[0].Y + Phi2*self.Points[1].Y + Phi3*self.Points[2].Y
			Z=Phi1*self.Points[0].Z + Phi2*self.Points[1].Z + Phi3*self.Points[2].Z
			Param[2].append([X,Y,Z])

			Param[0].append(1./6)
			Param[1].append((1./6, 4./6))
			Phi1=Phi_chap(1, 1./6, 4./6)
			Phi2=Phi_chap(2, 1./6, 4./6)
			Phi3=Phi_chap(3, 1./6, 4./6)
			X=Phi1*self.Points[0].X + Phi2*self.Points[1].X + Phi3*self.Points[2].X
			Y=Phi1*self.Points[0].Y + Phi2*self.Points[1].Y + Phi3*self.Points[2].Y
			Z=Phi1*self.Points[0].Z + Phi2*self.Points[1].Z + Phi3*self.Points[2].Z
			Param[2].append([X,Y,Z])

		return Param




def Phi_chap( i, e, n):
		if i==1:
			return 1 - e - n
		elif i==2:
			return e
		else:
			return n

def phiRef( i, param):
	Res=[]
	for j in range(len(param)):
		Res.append(Phi_chap(i, param[j][0], param[j][1]))
	return Res



class Mesh:
	def __init__(self):
		

		"""for p in Points:
			self.Points.append(copy.deepcopy(p))

		for S in Segments:
			self.Segments.append(copy.deepcopy(S))

		for T in Triangles:
			self.Triangles.append(copy.deepcopy(T))"""


	def GmshToMesh(self, filename=None):
		gmsh.initialize(sys.argv)
		if filename:
			Mesh=gmsh.merge(filename)

		(nodeTag, coord, paramcoord)=gmsh.model.mesh.getNodes() #nodeTag: Liste des tags des points
																#coord: Liste des coord des points de la forme [xp1, yp1, zp1, xp2, yp2...]
		self.Points=[0 for x in range(len(nodeTag))]
		self.Segments=[]
		self.Triangles=[]
		self.Tag1DId=[]														
		self.Tag2DId=[]
		Tag1DIdlen=0 # [(Phyical_Tag1 de dim 1, [indice_Segment1, ..., indice_Segmentn]), ..., (Physical_TagT, [indice_Segment1, ..., indice_Segmentm])]
		Tag2DIdlen=0# [(Phyical_Tag1 de dim 2, [indice_Triangle1, ..., indice_Trianglen]), ..., (Physical_TagT, [indice_Triangle1, ..., indice_Trianglem])]
		#print(nodeTag)
		for i in range(len(nodeTag)):
			#print(nodeTag[i])
			self.Points[nodeTag[i] -1]=Point(nodeTag[i], coord[3*i], coord[3*i+1], coord[3*i+2]) #Les points sont à l'indice de leur physical_tag -1

		#print(len(self.Points))

		dimTags=gmsh.model.getPhysicalGroups() #List de paires (dim, tag)
		LendimTags=len(dimTags)
		#print(dimTags)
		for i in range(LendimTags):
			dim=dimTags[i][0]

			if dim==1:
				self.Tag1DId.append((dimTags[i][1], []))
				Tag1DIdlen +=1
			if dim == 2:
				self.Tag2DId.append((dimTags[i][1], []))
				Tag2DIdlen +=1

			Tags=gmsh.model.getEntitiesForPhysicalGroup(dim, dimTags[i][1])
			#print(Tags)
			for j in range(len(Tags)):
				elementTypes, elementTags, nodeTags = gmsh.model.mesh.getElements(dim, Tags[j]) #elementTypes: liste des type de dimenssion: 1 si segment, 2 si triangles
																								#elementTags: Liste de liste: chaque sous liste contient les tags des sous-éléments constituant l'élément
																								#nodeTags: Liste de liste contenant les tags des points correspondant aux sous-éléments
				#print(elementTypes)
				for k in range(len(elementTypes)):
					if elementTypes[k]==1:
						for z in range(0,len(elementTags[k])):
							S=Segment([self.Points[nodeTags[k][2*z] -1], self.Points[nodeTags[k][2*z+1] -1]], elementTags[k][z])
							self.Segments.append(S)
							self.Tag1DId[Tag1DIdlen-1][1].append(S.id)

					if elementTypes[k]==2:
						for z in range(0,len(elementTags[k])):
							T=Triangle([self.Points[nodeTags[k][3*z] -1], self.Points[nodeTags[k][3*z+1] -1], self.Points[nodeTags[k][3*z+2] -1]], elementTags[k][z])
							self.Triangles.append(T)
							self.Tag2DId[Tag2DIdlen-1][1].append(T.id)



	def getElements(self,dim, physical_tag):
		Res=[]
		if dim==1:
			for S in self.Tag1DId:
				if S[0]==physical_tag:
					for x in S[1]:
						Res.append(self.Segments[x])
					break

		elif dim==2:
			for T in self.Tag2DId:
				if T[0]==physical_tag:
					for x in T[1]:
						#print(x)
						Res.append(self.Triangles[x])
					break

		return Res

	def getPoints(self, dim, physical_tag):
		Res=[]
		Lelement=self.getElements(dim, physical_tag)
		#print(len(Lelement))
		ListTag=[]
		for E in Lelement:
			for p in E.Points:
				if p.ID_Glob not in ListTag:
					Res.append(p)
					ListTag.append(p.ID_Glob)
		#print(ListTag)		
		return Res




"""
P1=Point(1,1,1,1)
P2=Point(2,0,1,0)
P3=Point(3,3,7,0)
S1=Segment([P1, P2], 1)
#print(S1.jac())
T1=Triangle([P1, P2, P3], 1)
#print(T1.jac())
#print(Point.identifiant)

M=Mesh([P1,P2,P3],[S1],[T1])
M=Mesh()
M.GmshToMesh("Fonction_Forme.msh")
print(M.getPoints(1,1)[0])
print(M.getPoints(1,1)[1])
#print(M.getElements(2,50)[0].Points[0])
print(M.Triangles[0].Points[0])
param=M.Triangles[0].gaussPoint()
print(Integrale(M,2,F))"""

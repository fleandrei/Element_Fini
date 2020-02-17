# Element_Fini

Projet en groupe: Guillaume Nguyen et Andrei Fleiser

Ce projet contient les fichiers code suivants:

	- "Maillage.py": Contient les classes et méthodes relatives à la gestion du maillage (classes Points, Triangles, Mesh...). On y trouve notament la fonction "gener_Gmsh_carre" qui génére la 		géométrie et le maillage gmsh (Question 1) nécéssaire à la suite
	
	- "Assemblage_Matrice.py": Contient toutes les fonctions et classes relatives à l'utilisatio de l'algèbre linéaires (classe Triplet, fonctions calculant les matrices de masse, de rigidité , le 	vecteur B...)

	- "Function.py": Contient diverses fonctions utiles dont notamment la fonction permettant de prendre en compte la condition de Dirichlet




Ce projet contient également les script exécutables suivants:
	
	- "Main1.py": Script répondant à la question 2. Il se lance en ligne de commande avec l'interpréteur python 2.7.15 (Le code ne compile pas en python 3), sans aucun argument et calcule la solution 		pour h=0.25. Le programme affiche le graphique de la solution approchée, ainsi que la solution exacte (l'interpolation linéaire). L'image de la solution approchée est sauvegardée dans l'image 	"Solution_Approchée"

	- "Main2.py": Script répondant à la question 3 et qui calcule l'erreur logarithmique en norme L2 entre la solution exacte et la solution approchée. Il existe 2 façons de lancer ce script:

		- Avec 1 argument: "h": Il s'agit de la valeur de h. Le script va alors afficher l'image de la solution approchée correspondante. Il affichera également dans le terminal l'erreur trouvée.
	
		- Avec 3 arguments: "Hmin Hmax pas" : Cela permet de calculer l'erreur pour plusieurs valeurs de h comprises entre "Hmin" et "Hmax" et espacées par un interval de "pas". Le script renvoie 			l'image de la courbe correspondante "Evolution_Erreur"; ainsi qu'un fichier Json "Résultat_Erreur" avec les données correspondantes 



Notes: 

La solution approchée obtenue suit assez bien les variations de la solution exacte. Les zones où la fonction est maximale coincident et il en est de même pour les zones où la fonction est minimale. 
Malheureusement, l'amplitudes de la fonction approchée est plus grande d'un facteur environ égal à (1+2*Pi*Pi). On réalise ainsi que si la solution exacte avait été la fonction 'f', notre solution approchée aurrait été adéquate. 
Par ailleur, étant donné qu'on calcule l'erreur comme étant la norme L2 de la différence des vecteurs U et Uref (vecteur contenant les valeurs de la fonction exacte pour les sommets du maillage), plus le nombre de points du maillage augmente et plus l'erreur devient importante (car on somme plus de différences de nombres ayant des normes très différentes). Ainsi contrairement à ce qui devrait se passer, l'erreur augmente quand h diminue. Dans ces circompstances, calculer la vitesse de convergence n'a pas de sens.
La grande norme des valeurs de U semblent venir du fait que les valeurs de certains coefficients de l'inverse de A sont eux mêmes assez élevées. Pourtant, les matrices de masse et de rigidité ont bien passé les testes de validités proposés lors du TP (Vtransposé*M*V=1,   K*V=0 (vecteur nulle) ; pour V=(1,...,1) ). 

A titre purement indicatif, nous avons voulu voir ce que cela aurrait donné si la fonction 'f' avait été la solution exacte. Vous pouvez trouver ces résultats dans l'image "Evolution_Erreur_Sol_F" et dans le fichier au format Json "Résultat_Erreur_Sol_F"

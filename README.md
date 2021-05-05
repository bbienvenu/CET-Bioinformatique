# CET-Bioinformatique - Projet individuel

  
<p align="justify">
  
  ## Algorithme de Needleman et Wunsch : détermination de tous les alignements optimaux
  
</p>

L’objectif du projet est d’écrire un programme implémentant l’algorithme de Needleman et Wunsch permettant de déterminer **tous** les alignements optimaux.

On essaiera de suivre la démarche suivante :

- construction du graphe des chemins optimaux : 

Chaque case de la matrice peut être décrite par un couple d’entiers (i,j) où i et j sont les indices de la ligne et de la colonne de la matrice contenant cette case.
Si m et n désignent la longueur des deux séquences, la matrice a une taille (m+1)*(n+1). En allant de la case d’indice (m,n) à la case d’indice (0,0), on peut construire un graphe 
des cases décrivant l’ensemble des chemins optimaux. Chaque case sera représentée par un nœud ayant pour fils le ou les nœud(s) représentant la ou les case(s) voisine(s) 
ayant permis d’obtenir le score optimal de cette case. Le graphe obtenu pour l’exemple précédent est présenté ci-contre. Ce graphe peut être construit récursivement en partant 
de la case d’indice (m,n).

- détermination de tous les chemins optimaux :

En parcourant le graphe des chemins optimaux, chaque chemin optimal peut être déterminé. La liste des chemins optimaux dans l’exemple précédent est la suivante :

(0, 0), (1, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 6), (7, 6), (8, 7), (9, 8)
(0, 0), (1, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 6), (7, 7), (8, 7), (9, 8)
(0, 0), (1, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 7), (8, 7), (9, 8)
(0, 0), (1, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 6), (7, 7), (8, 8), (9, 8)
(0, 0), (1, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 7), (8, 8), (9, 8)
(0, 0), (1, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 8), (9, 8)

Cette liste peut être construite récursivement à partir du graphe.

- détermination de tous les alignements optimaux : 

Le parcours de chacun des chemins optimaux permet de construire l’alignement optimal correspondant.

# -*- coding: utf-8 -*-

#-------------------------------------------------------------------- #
# algorithme de Needleman et Wunsch
#
# cours : CET bioinformatique
# Projet : détermination de tous les alignements optimaux issus de l’algorithme de Needleman et Wunsch
# auteur : Bienvenu Bambi
#-------------------------------------------------------------------- #


# ---------------  calcul de la matrice de distances  --------------- #


def match_ou_mismatch(res1, res2, match, mismatch):
    '''
    fonction renvoyant le score d'alignement de deux résidus
    '''
    
    if res1 == res2:
        return match
    else:
        return mismatch


def calcul_scores(seq1, seq2, match, mismatch, gap):
    '''
    fonction renvoyant la matrice des scores d'alignement optimal de
    tous les préfixes d'une séquence avec tous les préfixes d'une autre
    séquence
    '''
    
    # longueur des séquences
    lg_seq1 = len(seq1)
    lg_seq2 = len(seq2)
    
    # création et initialisation de la matrice de scores
    scores = list()
    for i in range(lg_seq1+1):
        scores.append([0]*(lg_seq2+1))
    
    # calcul des scores de la première colonne
    for i in range(1,lg_seq1+1):
        scores[i][0] = i * gap
    
    # calcul des scores de la première ligne
    for j in range(1,lg_seq2+1):
        scores[0][j] = j * gap
        
    # calcul des autres scores
    for i in range(1, lg_seq1+1):
        for j in range(1, lg_seq2+1):
            score_gauche = scores[i][j-1] + gap
            score_diagonale = scores[i-1][j-1] + match_ou_mismatch(seq1[i-1], 
                                                                   seq2[j-1],
                                                                   match,
                                                                   mismatch)
            score_haut = scores[i-1][j] + gap
            scores[i][j] = min(score_gauche, score_diagonale, score_haut)
    
    # renvoi de la matrice
    return scores



# --------------- affichage de la matrice de distances -------------- #


def affichage_scores(seq1, seq2, scores):
    '''
    fonction renvoyant une chaîne de caractères permettant l'affichage 
    d'une matrice de scores
    '''
    
    # initialisation de la chaine de caractères contenant le resultat
    res = ''

    # écriture de la séquence 2
    res += ' ' * 2
    res += ' ' * 4
    for c in seq2:
        res += ' ' * 3
        res += c
    res += '\n'
    
    # écriture de la séquence 1 et de la matrice
    for i in range(len(seq1)+1):
        # écriture du résidu de la séquence 1
        if i > 0:
            res += seq1[i-1]
            res += ' '
        else:
            res +=  ' ' * 2
        # écriture de la ligne de la matrice
        for j in range(len(seq2)+1):
            res += '%4s' % (scores[i][j],)
        res += '\n'
    
    # renvoi de la chaine de caractères contenant le résultat
    return res



# ---------------  recherche d'un alignement optimal  --------------- #


def calcul_un_alignement(seq1, seq2, match, mismatch, gap, scores):
    '''
    fonction renvoyant l'un des alignements optimaux obtenus à partir 
    d'une matrice de scores
    '''
    
    # séquences alignées
    seq1_aln = ''
    seq2_aln = ''
    
    # on commence l'alignement par la fin
    i = len(seq1)
    j = len(seq2)
    
    # parcours de la matrice
    while i > 0 and j > 0:
        # gap dans la séquence 1
        if scores[i][j] == scores[i][j-1] + gap:
            seq1_aln = '-' + seq1_aln
            seq2_aln = seq2[j-1] + seq2_aln
            j -= 1
        # match ou mismatch
        elif scores[i][j] == scores[i-1][j-1] + match_ou_mismatch(seq1[i-1],
                                                                  seq2[j-1],
                                                                  match,
                                                                  mismatch):
            seq1_aln = seq1[i-1] + seq1_aln
            seq2_aln = seq2[j-1] + seq2_aln
            i -= 1
            j -= 1
        # gap dans la séquence 2
        else:
            seq1_aln = seq1[i-1] + seq1_aln
            seq2_aln = '-' + seq2_aln
            i -= 1
    
    # s'il y a un gap au début de la deuxième séquence
    while i > 0:
        seq1_aln = seq1[i-1] + seq1_aln
        seq2_aln = '-' + seq2_aln
        i -= 1
    
    # s'il y a un gap au début de la première séquence
    while j > 0:
        seq1_aln = '-' + seq1_aln
        seq2_aln = seq2[j-1] + seq2_aln
        j -= 1
    
    # renvoi des séquences alignées
    return seq1_aln, seq2_aln



# ---------------  recherche de tous les alignements optimaux  --------------- #


# construction du graphe des chemins optimaux

# classe représentant le graphe
class Tree:
    def __init__(self, key):
        self.gauche = None
        self.diagonale = None
        self.haut = None
        self.value = key

    def __str__(self):
        node = self
        if node != None:
            s = str(node.value)
            if node.gauche: 
                s = node.gauche.__str__() + s
            if node.diagonale:
                s += ' ' + node.diagonale.__str__()
            if node.haut:
                s += ' ' + node.haut.__str__()
            return s
        return ''


def calcul_alignements(seq1, seq2, match, mismatch, gap, scores, i=None, j=None, graphe=None):
    '''
    Fonction renvoyant le graphe des chemins optimaux.
    Chaque case de la matrice des scores est représentée  
    par un nœud ayant pour fils le ou les nœud(s)représentant 
    la ou les case(s) voisine(s) ayant permis d’obtenir le 
    score optimal de cette case.
    '''

    if (i==0 and j==0):
        return

    if i == 0 and j != 0:
        graphe.gauche = Tree((i,j-1))
        calcul_alignements(seq1, seq2, match, mismatch, gap, scores, i, j-1, graphe.gauche)

    if j == 0 and i != 0:
        graphe.haut = Tree((i-1,j))
        calcul_alignements(seq1, seq2, match, mismatch, gap, scores, i-1, j, graphe.haut)

    # Initialisation
    if (i==None and j==None and graphe==None):
        i = len(seq1)
        j = len(seq2)
        graphe = Tree((i,j))

    # gap dans la séquence 1
    if scores[i][j] == scores[i][j-1] + gap:
        graphe.gauche = Tree((i,j-1))
        calcul_alignements(seq1, seq2, match, mismatch, gap, scores, i, j-1, graphe.gauche)

    # match ou mismatch
    if scores[i][j] == scores[i-1][j-1] + match_ou_mismatch(seq1[i-1],
                                                                seq2[j-1],
                                                                match,
                                                                mismatch):
        graphe.diagonale = Tree((i-1,j-1))
        calcul_alignements(seq1, seq2, match, mismatch, gap, scores, i-1, j-1, graphe.diagonale)

    # gap dans la séquence 2
    if scores[i][j] == scores[i-1][j] + gap:
        graphe.haut = Tree((i-1,j))
        calcul_alignements(seq1, seq2, match, mismatch, gap, scores, i-1, j, graphe.haut)

    return graphe


# détermination de tous les chemins optimaux

def find_optimal_paths(tree):
    '''
    Fonction renvoyant la liste des chemins optimaux 
    obtenus à partir du graphe préalablement construit.
    '''

    # Initialisation de la liste des chemins
    paths = list()

    if not (tree.gauche or tree.diagonale or tree.haut):
        return [[tree.value]]

    if tree.gauche:
        paths.extend([[tree.value] + child for child in find_optimal_paths(tree.gauche)])

    if tree.diagonale:
        paths.extend([[tree.value] + child for child in find_optimal_paths(tree.diagonale)])

    if tree.haut:
        paths.extend([[tree.value] + child for child in find_optimal_paths(tree.haut)])

    return paths


# détermination de tous les alignements optimaux

def find_optimal_alignments(seq1, seq2, paths):
    '''
    Fonction renvoyant tous les alignements optimaux 
    en parcourant chacun des chemins optimaux trouvés avant.
    '''

    # Initialisation de la liste des alignements optimaux
    alignments = list()

    # on remet la liste dans le bon ordre (on commence à partir de l'indice (0,0))
    for liste in paths:
        liste.reverse()
    
    for path in paths:
        seq1_aln = str()
        seq2_aln = str()
        n = len(path)

        for i in range(1,n): # on commence le range à 1 puisque (0,0) sera toujours le premier (dernier) élément
            if (path[i-1][0] == path[i][0]-1) and (path[i-1][1] == path[i][1]-1):
                seq1_aln += seq1[path[i][0]-1]
                seq2_aln += seq2[path[i][1]-1]

            elif (path[i-1][0] == path[i][0]-1) and (path[i-1][1] == path[i][1]):
                seq1_aln += seq1[path[i][0]-1]
                seq2_aln += '-'

            elif (path[i-1][0] == path[i][0]) and (path[i-1][1] == path[i][1]-1):
                seq1_aln += '-'
                seq2_aln += seq2[path[i][1]-1]

        alignments.append((seq1_aln, seq2_aln))

    return alignments



# ---------------         programme principal         --------------- #


    # --- données --- #

# séquences
data = list()
data.append(('ATGCATTTT', 'ACTGCATT'))       # exemple de l'étape 1
data.append(('ATGC', 'ATCCGC'))              # exemple 1 (= exemple du cours)
data.append(('CGATGGC', 'TACGATGGC'))        # exemple 2
data.append(('CGATCGTGACTT', 'ATCGTGAC'))    # exemple 3

# choix des séquences
sequences = data[0]
sequence_1 = sequences[0]
sequence_2 = sequences[1]

# paramètres de l'algorithme
match = 0
mismatch = 1
gap = 2


# --- calculs --- #

# calcul de la matrice de scores
scores = calcul_scores(sequence_1, 
                       sequence_2,
                       match,
                       mismatch,
                       gap)

# construction du graphe des chemins optimaux
tree = calcul_alignements(sequence_1, sequence_2, match, mismatch, gap, scores, i=None, j=None, graphe=None)

# détermination de tous les chemins optimaux
paths = find_optimal_paths(tree)

# détermination de tous les alignements optimaux
alignements = find_optimal_alignments(sequence_1, sequence_2, paths)

# --- affichage du résultat --- #

# séparateur
separateur = '_'*30 + '\n'

# titre
print('')
print('*'*52)
print('Alignement selon l\'algorithme de Needleman et Wunsch')
print('*'*52)

print('%s' % (separateur,))

# paramètres de l'algorithme
print('Paramètres de l\'algorithme :')
print('- %-25s%2s' % ('score d\'un match : ', match))
print('- %-25s%2s' % ('score d\'un mismatch : ', mismatch))
print('- %-25s%2s' % ('score d\'un gap : ', gap,))

print('%s' % (separateur,))

# séquences à aligner
print('Séquences à aligner :')
print('séquence 1 : longueur = %s' % (len(sequence_1),))
print('%s' % (sequence_1,))
print('séquence 2 : longueur = %s' % (len(sequence_2),))
print('%s' % (sequence_2,))

print('%s' % (separateur,))

# matrice de scores
print('Matrice de scores :')
print('%s' % (affichage_scores(sequence_1, sequence_2, scores)[:-1]))

print('%s' % (separateur,))

# alignement

n = len(alignements)
for i in range(n):
    print('Alignement %s :' % str(i+1))
    print(alignements[i][0])
    print(alignements[i][1])
    print('%s' % (separateur,))

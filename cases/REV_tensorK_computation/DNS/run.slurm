#!/bin/bash

#SBATCH -N 1
# nombre du noeuds

#SBATCH -n 20
# nombre de taches MPI

#SBATCH --ntasks-per-node=20
# Nombre de tâches MPI par noeud (doit correspondre au nombre de total
# de taches MPI divisé par le nombre de noeuds)

#SBATCH --ntasks-per-core=1
# nimbre de taches MPI per noeuds (valeur a laisser a 1)


# Mémoire réservée par noeud en MB (Si cette option n’est pas spécifiée
# en en-tête du script, la totalité de la mémoire du noeud est allouée.)

#SBATCH -o out.log
# Le fichier dans lequel sera écrit les sorties

#SBATCH -e err.log
# Le fichier dans lequel sera écrit les messages d’erreur

#SBATCH --time=20:00:00 
# Le temps maximum pour lequel le job va tourner (Passé ce délai, le job est automatiquement arrêté.)

#SBATCH --mail-user=nluminar@imft.fr
# Adresse mail à laquelle recevoir les notifications

#SBATCH --mail-type=ALL
# Cas de figure dans lesquels le mail est envoyé : BEGIN, END, FAIL, REQUEUE,
# and ALL. (Si cette option n’est pas spécifiée en en-tête du script,
# aucun mail ne sera envoyé.)

#-------------------------------------------------------

#pyFoamPrepareCase.py . --parameter-file=parameters.basic 2>&1 > prepare.log
#pyFoamPrepareCase.py . --no-mesh-create  --parameter-file=parameters.basic
mpirun -n 20 icoGFoam -parallel

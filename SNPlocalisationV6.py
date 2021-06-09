# coding=utf-8


__author__ = 'ericfournier'

"""
Exemple de commande dans /home/eric/ProjetPersonnel/Programmation/ProjetPython/NGS/NotForGalaxy
python SNPlocalisationV6.py -g ref.gb -s prefix-profile_test.tsv  -o MyRes2.txt --stat MyStat.txt

Modifications par rapport a SNPlocalisationV2 :
    - les positions des snv sont trie par ordre croissant dans le fichier de sortie
    - statistiques sommaires : nombre de snv dans chaque feature + longueur totale de chaque feature dans la reference

Note: un element sur le brin (+) peut chevaucher un autre element sur le brin (-)




"""


import re
import os
import sys
from Bio import SeqIO
import logging
import argparse
import threading
import time

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()
logger.handlers = []  # This is the key thing for the question!
handler = logging.StreamHandler()
formatter = logging.Formatter("%(asctime)s: %(message)s")
handler.terminator = ""
handler.setFormatter(formatter)
logger.addHandler(handler)


global running
running = True

def thread_function(name):
    """
    Message d attente
    :param name:
    :return:
    """

    while running:
        time.sleep(1)
        logging.info("Calcul en cours")
    logging.info("Terminé")


x = threading.Thread(target=thread_function, args=(1,))
x.start()


def GetIntergenicLength(locations,total_seq_length):
    """
    Calculer la longueur total des regions intergenic
    :param locations:
    :param total_seq_length:
    :return:
    """

    #region intergenic en 5'
    start_intergenic_length = locations[0].start - 1

    #region intergenic en 3'
    end_intergenic_length = total_seq_length - locations[len(locations)-1].end

    total_intergenic = start_intergenic_length + end_intergenic_length
    index = 0
    for location in locations[:-1]:
        if(location.end < locations[index+1].start):
            total_intergenic += locations[index+1].start - location.end - 1
        index += 1

    return total_intergenic

def ComputeStat(data,features):
    """
    Compter le nombre de snv dans chacun des features
    :param data:
    :param features:
    :return:
    """
    global  feature_count
    feature_count = {}

    for feat in features:
        feature_count[feat] = 0

    for feat in dict(data).values():
        feature_count[feat[0]] += 1

def ExtractReferenceFeature(gb_file):
    """
    Obtenir la liste des features dans la reference et
    calculer la longueur de chacun des features dans la reference
    :param gb_file:
    :return:
    """

    features = set([])
    feature_location = []
    global feature_length
    feature_length = {}
    for feat in reference_record.features:
        if feat.type not in ['gene','source']:

            features.add(feat.type)

            if feat.type not in feature_length.keys():
                feature_length[feat.type] = len(feat.location)
            else:
                feature_length[feat.type] += len(feat.location)
            feature_location.append(feat.location)

    feature_length['intergenic'] = GetIntergenicLength(feature_location,len(gb_file.seq))

    features.add('intergenic')
    return features

#Les options de la ligne de commande
parser = argparse.ArgumentParser(description="Localisateur de SNV")
parser.add_argument("-g","--gbfile", help="Path vers le fichier genbank de référence",required=True)
parser.add_argument("-s", "--snvfile",help="Path vers le fichier de SNV",required=True)
parser.add_argument("-o", "--out",help="Path vers le fichier de sortie",required=True)
parser.add_argument("--stat",help="Path vers le fichier statistique",required=True)

args_commandline = parser.parse_args(sys.argv[1:])
args = args_commandline.__dict__
gb_file = args["gbfile"]
snv_file = args["snvfile"]
out_file = args["out"]
stat_file = args["stat"]

#On s assure que les path existent
if not (os.path.exists(gb_file)) or not(os.path.exists(snv_file)):
    logging.error("L'un des fichiers input est inexistant")
    exit(1)

if str(out_file).count('/') > 0:
    if not(os.path.exists(os.path.dirname(out_file))):
        logging.error("Le répertoire vers le fichier de sortie est inexistant")
        exit(1)

if str(stat_file).count('/') > 0:
    if not(os.path.exists(os.path.dirname(stat_file))):
        logging.error("Le répertoire vers le fichier statistique est inexistant")
        exit(1)

#Le fichier de reference
reference_record = SeqIO.read(gb_file,'genbank')

#Liste des features dans la reference
feature_set = ExtractReferenceFeature(reference_record)

pos_set=set() # set de snv du fichier de SNV
pos_list=[] # liste de snv du fichier d SNV
pos_dict={} # dictionnaire de type key=position snv value = [feature,product]

with open(snv_file) as readf:  # on extraite les snv positions du fichier -profil.tsv
    pos_list = readf.readline().split('\t')[1:]
    pos_list = [re.search(r'^.+_(\d+)', gi).group(1) for gi in pos_list]
readf.close()

for pos in pos_list:
    pos = int(pos)
    pos_dict[pos] = []
    pos_set.add(pos)

del(pos_list)

#On construit le pos_dict
for feat in reference_record.features:
    if feat.type in feature_set:
        for pos in pos_set:
            if int(pos) in feat and len(pos_dict[pos]) == 0:
                pos = int(pos)
                pos_dict[pos].append(feat.type)
                try:
                    pos_dict[pos].append(feat.qualifiers.get('product'))
                except:
                    pos_dict[pos].append('NA')


print pos_dict

#snv intergenic
for pos in pos_set:
    if len(pos_dict[pos]) == 0:
        pos_dict[pos].append('intergenic')
        pos_dict[pos].append('NA')

del(pos_set)

#fichier contenant les positions des snv
out_handle=open(out_file,'w')
out_handle.write('Gene\tFeature\tPosition\n') # le header
for pos in sorted(pos_dict.keys()):  # on ecrit le product tab le feature tab la position du snp
    out_handle.write(str(pos_dict[pos][1]) + '\t' + str(pos_dict[pos][0]) + '\t' + str(pos) + '\n')
out_handle.close()

#Calcul des statistiques sommaires
ComputeStat(pos_dict,feature_set)
stat_handle = open(stat_file,'w')
stat_handle.write('SNV in Regions\n')
for feat in sorted(feature_count.keys()):
    stat_handle.write('\t' + feat + " : " +str(feature_count[feat]) + '\n')
stat_handle.write('\nNucleotides in Regions\n')
for feat in sorted(feature_length.keys()):
    stat_handle.write('\t' + feat + " : "+ str(feature_length[feat]) + '\n')

stat_handle.close()
running = False






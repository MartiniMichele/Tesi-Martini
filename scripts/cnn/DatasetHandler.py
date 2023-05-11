import math
import os
from pathlib import Path
import shutil
import random
import glob

source_path = Path(__file__).resolve()
source_dir = Path(source_path.parent.parent.parent)

database = "lsu"
livello = "eukaryota"
completo = "%s_%s" % (database, livello)
sottoclasse = "eukaryota"

dataset_dir = Path(str(source_dir) + "/Classification/%s_DATASET/" % database.upper())
data_dir = Path(str(source_dir) + "/Classification/IMMAGINI DA DIVIDERE/%s/" % (completo))
train_dir = Path(str(source_dir) + "/Classification/IMMAGINI DA DIVIDERE/%s/%s/train/" % (completo, sottoclasse))
valid_dir = Path(str(source_dir) + "/Classification/IMMAGINI DA DIVIDERE/%s/%s/valid/" % (completo, sottoclasse))
test_dir = Path(str(source_dir) + "/Classification/IMMAGINI DA DIVIDERE/%s/%s/test/" % (completo, sottoclasse))

if os.path.isdir(train_dir) is False: os.makedirs(train_dir)
if os.path.isdir(valid_dir) is False: os.makedirs(valid_dir)
if os.path.isdir(test_dir) is False: os.makedirs(test_dir)

filelist = os.listdir(data_dir)
sub_dir_path = []
'''
cerca tutte le sottocartelle relative al phylum(o superkingdom) e poi itera al loro interno per spostare i file secondo
la proporzione 70/20/10 nelle rispettive cartelle train/valid/test
'''
if len(os.listdir(train_dir)) == 0:
    for x in filelist:
        if x != "DATASET":
            sub_dir_path.append(str(data_dir) + '/' + x)
    for x in sub_dir_path:
        print(str(x))

    for i in sub_dir_path:
        filelist2 = os.listdir(i)
        train_part = math.floor(len(filelist2) * 0.7)
        valid_part = math.floor(len(filelist2) * 0.2)
        test_part = math.floor(len(filelist2) * 0.1) - 1
        print("phylum: " + str(i))
        print("70: " + str(train_part))
        print("20: " + str(valid_part))
        print("10: " + str(test_part))

        os.chdir(i)
        for j in random.sample(glob.glob('*.png'), train_part):
            source_path = i + '/' + j
            shutil.move(source_path, train_dir)
        for j in random.sample(glob.glob('*.png'), valid_part):
            source_path = i + '/' + j
            shutil.move(source_path, valid_dir)
        for j in random.sample(glob.glob('*.png'), test_part):
            source_path = i + '/' + j
            shutil.move(source_path, test_dir)


'''
controlla se è presente la cartella "Actinobacteria"(può variare) comune a tutte e 3 le cartelle(train, valid, test) e quindi sicuramente presente
se si va avanti e crea le sottocartelle(labels) con i nomi dei phylum presenti.
Per creare le cartelle dei labels vengono letti tutti i file .png, viene estratta la prima parte del nome relativa al phylum
e poi ci viene spostato il file
'''
if(os.path.isdir(dataset_dir) is False):
    os.makedirs(dataset_dir)
os.chdir(train_dir)
for file in os.listdir():
    if file.endswith('.png'):
        phylum = file.split('_')[0]
        if os.path.isdir(phylum) is False:
            os.makedirs(phylum)
            print("CARTELLA TRAIN CREATA")
        dir_path = str(dataset_dir) + '/train/' + phylum
        print(f"il file è {str(file)}")
        if(os.path.isdir(dir_path) is False):
            os.makedirs(dir_path)
            print("dir_path creato")
        shutil.move(file, dir_path)

os.chdir(valid_dir)
for file in os.listdir():
    if file.endswith('.png'):
        phylum = file.split('_')[0]
        if os.path.isdir(phylum) is False:
            os.makedirs(phylum)
            print("CARTELLA VALID CREATA")
        dir_path = str(dataset_dir) + '/valid/' + phylum
        print(f"il file è {str(file)}")
        if (os.path.isdir(dir_path) is False):
            os.makedirs(dir_path)
            print("dir_path creato")
        shutil.move(file, dir_path)

os.chdir(test_dir)
for file in os.listdir():
    if file.endswith('.png'):
        phylum = file.split('_')[0]
        if os.path.isdir(phylum) is False:
            os.makedirs(phylum)
            print("CARTELLA TEST CREATA")
        dir_path = str(dataset_dir) + '/test/' + phylum
        print(f"il file è {str(file)}")
        if (os.path.isdir(dir_path) is False):
            os.makedirs(dir_path)
            print("dir_path creato")
        shutil.move(file, dir_path)

from scripts.cgr import CGRepresentation, FCGR
import os
from pathlib import Path
from Bio.Seq import MutableSeq
import pandas as pd

from scripts.cgr.FCGR import FrequencyCGR
from Levenshtein import distance as levenshtein_distance

'''
Questa classe ha il compito di gestire la generazione di immagini CGR tramite l'utilizzo del codice della libreria
adkwazar/CGR
'''


class CGRHandler:
    source_path = Path(__file__).resolve()
    source_dir = Path(source_path.parent.parent.parent)

    def __init__(self, CGR_type, outer_representation, rna_2structure, data_dir, save_dir):
        self.CGR_type = CGR_type
        self.outer_representation = outer_representation
        self.rna_2structure = rna_2structure
        self.sequence = None
        self.source = data_dir
        self.data_dir = data_dir
        self.save_dir = save_dir

    '''
    Inizializza le cartelle utili alla CNN
    '''

    def init_dirs(self, save_dir):
        # cartelle urilizzate per l'esperimento
        self.data_dir = Path(str(self.source_dir) + "/FASTA/%s/" % self.data_dir)
        cgr_dir = Path(str(self.source_dir) + save_dir)

        if os.path.isdir(cgr_dir) is False:
            os.makedirs(cgr_dir)
        self.save_dir = Path(str(cgr_dir) + "/%s" % self.save_dir)

        os.chdir(cgr_dir)
        if os.path.isdir(self.save_dir) is False:
            os.makedirs(self.save_dir)

    '''
    Questo metodo scorre i file fasta nella cartella selezionata e ne estrae la sequenza
    '''

    def read_files(self, isFCGR, k):

        if isFCGR:
            self.init_dirs("/IMMAGINI_FCGR")
        else:
            self.init_dirs("/IMMAGINI_CGR")

        # Folder Path
        source_path = Path(__file__).resolve()
        source_dir = Path(source_path.parent.parent.parent)
        path = Path(str(source_dir) + "/FASTA/%s" % self.source)
        counter = 1
        # Change the directory
        os.chdir(path)

        # Read text File

        def read_fasta_file(file_path):
            with open(file_path, 'r') as f:
                line = f.readline()
                self.sequence = line.replace("\n", "")

        # iterate through all file
        for file in os.listdir():
            # Check whether file is in text format or not
            if file.endswith(".fasta"):
                file_path = f"{path}\{file}"

                # call read text file function
                read_fasta_file(file_path)
                print(file_path)
                self.generate_dataset(counter, isFCGR, k)
                counter += 1

    '''
    Questo metodo filtra le sequenze rimuovendo i caratteri che causano ambiguità
    '''

    def filter_sequence(self, sequence):

        filtered_sequence = MutableSeq(str(sequence))
        chars = ["Y", "N", "R", "M", "S", "W", "K", "D", "V", "B", "H", "P", "O"]

        if any(x in sequence for x in chars):

            char_count = self.count_char(sequence)

            for i in range(0, char_count):

                if filtered_sequence.find("Y") != -1:
                    filtered_sequence = filtered_sequence.replace("Y", "")

                elif filtered_sequence.find("N") != -1:
                    filtered_sequence = filtered_sequence.replace("N", "")

                elif filtered_sequence.find("R") != -1:
                    filtered_sequence = filtered_sequence.replace("R", "")

                elif filtered_sequence.find("M") != -1:
                    filtered_sequence = filtered_sequence.replace("M", "")

                elif filtered_sequence.find("S") != -1:
                    filtered_sequence = filtered_sequence.replace("S", "")

                elif filtered_sequence.find("W") != -1:
                    filtered_sequence = filtered_sequence.replace("W", "")

                elif filtered_sequence.find("D") != -1:
                    filtered_sequence = filtered_sequence.replace("D", "")

                elif filtered_sequence.find("B") != -1:
                    filtered_sequence = filtered_sequence.replace("B", "")

                elif filtered_sequence.find("H") != -1:
                    filtered_sequence = filtered_sequence.replace("H", "")

                elif filtered_sequence.find("V") != -1:
                    filtered_sequence = filtered_sequence.replace("V", "")

                elif filtered_sequence.find("K") != -1:
                    filtered_sequence = filtered_sequence.replace("K", "")

                elif filtered_sequence.find("P") != -1:
                    filtered_sequence = filtered_sequence.replace("P", "")

                elif filtered_sequence.find("O") != -1:
                    filtered_sequence = filtered_sequence.replace("O", "")

                else:
                    pass

        return filtered_sequence

    '''
    metodo di supporto a filter_sequence, ha il compito di contare il numero di caratteri ambigui presenti nella sequenza
    '''

    def count_char(self, sequence):
        tmp_sequence = MutableSeq(sequence)
        char_count = 0

        char_count += tmp_sequence.count("Y")
        char_count += tmp_sequence.count("N")
        char_count += tmp_sequence.count("R")
        char_count += tmp_sequence.count("M")
        char_count += tmp_sequence.count("S")
        char_count += tmp_sequence.count("W")
        char_count += tmp_sequence.count("K")
        char_count += tmp_sequence.count("D")
        char_count += tmp_sequence.count("H")
        char_count += tmp_sequence.count("V")
        char_count += tmp_sequence.count("B")
        char_count += tmp_sequence.count("P")
        char_count += tmp_sequence.count("O")

        return char_count

    '''Questo metodo ha il compito di generare le immagini CGR sfruttando la libreria adkwazar/CGR e filter_sequence 
    per filtrare le sequenze
    '''

    def generate_dataset(self, counter, isFCGR, k):

        if isFCGR:
            fcgr_sequence = self.filter_sequence(self.sequence)
            print("SEQUENZA UTILIZZATA: " + fcgr_sequence)
            print("COUNTER: " + str(counter))
            drawer = FrequencyCGR(fcgr_sequence)
            path = Path(str(self.save_dir) + "/FCGR_RNA_" + str(counter) + ".png")
            drawer.save_fcgr(k, path)
        else:
            bio_sequence = self.filter_sequence(self.sequence)
            print("SEQUENZA UTILIZZATA: " + bio_sequence)
            print("COUNTER: " + str(counter))
            drawer = CGRepresentation.CGR(bio_sequence, self.CGR_type, self.outer_representation, self.rna_2structure)
            drawer.representation()
            path = Path(str(self.save_dir) + "/CGR_RNA_" + str(counter) + ".png")
            drawer.plot(counter, path)

    def read_file(self, k_list, is_subsequence, csv_name=""):
        self.init_dirs("/IMMAGINI_FCGR")

        # Folder Path
        source_path = Path(__file__).resolve()
        source_dir = Path(source_path.parent.parent.parent)
        file_path = Path(str(source_dir) + "/FASTA/%s" % self.source)

        sequences = []

        with open(file_path, 'r') as f:
            tmp_str = f.name.split(".")[0].split("_")
            image_name = tmp_str[len(tmp_str) - 1]
            seq = ''
            for line in f:
                if line.startswith('>'):
                    if seq:
                        if len(sequences) < 10000:
                            sequences.append(seq)
                            seq = ''
                else:
                    seq += line.strip()
            if len(sequences) < 10000:
                sequences.append(seq)

        # no_duplicate_sequences = self.delete_duplicates(sequences)
        # print(f"Sono state trovate {len(no_duplicate_sequences)} sequenze diverse")
        self.generate_dataset_from_list(sequences, k_list, image_name, is_subsequence, csv_name)

    '''
    Metodo usato per comparare due sequenze
    '''

    def delete_duplicates_list(self, sequences):
        nd_sequences = [sequences[0]]

        for element in range(1, len(sequences)):
            for nd_element in nd_sequences:
                if not element == nd_element:
                    nd_sequences.append(element)
                else:
                    print(f"Sequenza {element} uguale a {nd_element}")
                    break

        return nd_sequences

    def delete_duplicates_string(self, nd_sequences, str_sequence):

        for element in nd_sequences:
            if element != str_sequence:
                nd_sequences.append(str_sequence)
                return str_sequence
            else:
                print(f"Sequenza {element} uguale a {str_sequence}")
                return False

    def generate_dataset_from_list(self, sequence_list, k_list, image_name, is_subsequence, csv_name=""):
        no_duplicates_list = []
        for i, sequence in enumerate(sequence_list):
            counter = i + 1
            print(f"L'elemento {counter} è: {sequence}")
            bio_sequence = self.filter_sequence(sequence)

            if is_subsequence:
                tmp_sequence = bio_sequence.replace("T", "U")
                correct_sequence = self.get_subsequence(i, tmp_sequence, csv_name)
            else:
                correct_sequence = bio_sequence.replace("T", "U")

            if i == 0:
                no_duplicates_list.append(correct_sequence)

            if self.delete_duplicates_string(no_duplicates_list, correct_sequence):
                print("SEQUENZA UTILIZZATA: " + correct_sequence)
                print("COUNTER: " + str(counter))

                drawer = FrequencyCGR(correct_sequence)
                # path = Path(str(self.save_dir) + f"/{image_name}_" + str(counter) + ".png")

                for k in k_list:

                    self.save_dir = self.save_dir if str(self.save_dir).endswith(str(k)) else Path(
                        str(self.save_dir).rsplit("_", 1)[0] + "_" + "K" + str(k))

                    if os.path.isdir(self.save_dir) is False:
                        os.makedirs(self.save_dir)
                    print("K-MER ATTUALE: " + str(k))
                    path = Path(str(self.save_dir) + f"/{image_name}_" + str(counter) + ".png")
                    drawer.save_fcgr(k, path)

    def get_subsequence(self, index, sequence, csv_name):
        sub_sequence = sequence[0:150]
        '''
        csv_path = Path(str(self.source_dir) + "/CSV/" + csv_name)
        file = pd.read_csv(csv_path, usecols=["molecola", "start", "finish"])
        sub_sequence = ""

        sequence_found = False
        for row in file.iloc:

            molecola = row.get("molecola")

            while sequence_found is False:

                if index == molecola:
                    start = row.get("start")
                    finish = row.get("finish")

                    sub_sequence = sequence[start - 1:finish]
                    sequence_found = True
            break
        '''
        return sub_sequence

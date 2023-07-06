from scripts.cgr.CGRHandler import CGRHandler
from scripts.cnn.CNN import CNN
from scripts.cnn.ResultPlotter import ResultPlotter
from pathlib import Path
import pandas as pd
import platform

'''
Metodo per l'implementazione della CLI, sfrutta altri metodi per semplificare la lettura
'''

file_list = ["RF00001_5S.fa", "RF00177_rRNA_SSU.fa", "RF02541_rRNA_LSU.fa"]
dataset_list = ["RF_2_DATASET_1K", "RF_2_DATASET_4K", "RF_2_DATASET_7K"]
kmer_list = [1, 4, 7]


def main():
    model_mk = 1
    batch_size = 32
    epochs = 30
    fl_filter = 32
    n_dropout = 1
    drop_value = 0.5
    n_layer = 3
    lr = 0.0001
    patience = 5

    print("SCEGLI IL PROCESSO DA ESEGUIRE:\n")
    print("1: GENERAZIONE IMMAGINI DA SEQUENZE RNA")
    print("2: CLASSIFICAZIONE IMMAGINI TRAMITE CNN")
    print("3: GENERAZIONE IMMAGINI E CLASSIFICAZIONE")
    process_input = input("la tua scelta: ")
    match process_input:
        case "1":
            imgen_case()
        case "2":
            cnn_case(model_mk, batch_size, epochs, fl_filter, n_dropout, drop_value, n_layer, lr, patience)
        case "3":
            auto_case(model_mk, batch_size, epochs, fl_filter, n_dropout, drop_value, n_layer, lr, patience)
            # imgen_case()
            # cnn_case(model_mk, batch_size, epochs, fl_filter, n_dropout, drop_value, n_layer, lr, patience)


'''
Metodo di supporto a main, contiene il codice per la scelta della generazione delle immagini
'''


def auto_case(model_mk, batch_size, epochs, fl_filter, n_dropout, drop_value, n_layer, lr, patience):
    '''
    for file in file_list:
        handler_istance = CGRHandler("RNA", False, False, file, file.split(".")[0] + "_images")
        handler_istance.read_file(kmer_list, True)

    #handler_istance = CGRHandler("RNA", False, False, "_LSU_rRNA_archaea - Copia.fa", "prova_image")
    #handler_istance.read_file(kmer_list, True, "sub_strings.csv")

    #csv_path = Path("C:/Users/Michele/Documents/GitHub/Tesi-Martini/CSV")
    #file = pd.read_csv(Path(str(csv_path) + "/sub_strings.csv"), usecols=["molecola", "start", "finish"])

    '''
    for dataset in dataset_list:
            dataset_directory = dataset.upper()
            k = int(dataset_directory.rsplit("_", 1)[1].replace("K", ""))
            print(f"DATASET UTILIZZATO: {dataset_directory}")

            cnn_instance = CNN(dataset_directory, model_mk, batch_size, epochs, fl_filter, 3, n_dropout, drop_value,
                               n_layer, lr, patience)

            print("RIEPILOGO RETE:")
            model = cnn_instance.create_model()
            print("----------GENERAZIONE BATCH DI DATI ED INIZIO TRAINING----------")
            datagen_list = cnn_instance.datagen()
            ########################
            ######## TRAIN
            ########################
            history = cnn_instance.train(model, datagen_list[0], datagen_list[1])
            save_train_val_to_txt(history, epochs, batch_size, lr, k, dataset)
            ########################
            ######## TEST
            ########################
            score = cnn_instance.test_evaluate(model, datagen_list[2])
            save_test_to_txt(score, epochs, batch_size, lr, k, dataset)


def imgen_case():
    fasta_directory = input(
        "\nINSERIRE IL NOME DELLA CARTELLA IN CUI SI TROVANO I FILE FASTA (sottocartella di "
        "'FASTA')").lower()
    isFCGR = input(
        "\nCHE TIPO DI IMMAGINI GENERARE? \n1: CGR \n2: FCGR \nla tua scelta: ")
    match isFCGR:
        case "1":
            images_directory = input(
                "\nINSERIRE IL NOME DELLA CARTELLA IN CUI SI SALVARE LE IMMAGINI (sottocartella di "
                "'IMMAGINI CGR')").lower()
            handler_istance = CGRHandler("RNA", False, False, fasta_directory, images_directory)
            handler_istance.read_files(False, 0)
        case "2":
            images_directory = input(
                "\nINSERIRE IL NOME DELLA CARTELLA IN CUI SI SALVARE LE IMMAGINI (sottocartella di "
                "'IMMAGINI FCGR')").lower()
            k = input("SCEGLIERE LUNGHEZZA DEI K-MER: ")

            is_subsequence = input("VUOI CONSIDERARE SOLAMENTE UNA SOTTOSEQUENZA? \nS: SI \n N: NO \nla tua scelta: ")

            match is_subsequence:
                case "S":
                    csv_name = input(
                        "\nINSERIRE IL NOME DEL FILE CSV (sottocartella di "
                        "'CSV')")
                    handler_istance = CGRHandler("RNA", False, False, fasta_directory, images_directory)
                    handler_istance.read_file(kmer_list, True, csv_name)
                case "N":
                    handler_istance = CGRHandler("RNA", False, False, fasta_directory, images_directory)
                    handler_istance.read_file(kmer_list, False)


'''
Metodo di supporto a main, contiene il codice per la scelta della classificazione delle immagini tramite CNN
'''


def cnn_case(model_mk, batch_size, epochs, fl_filter, n_dropout, drop_value, n_layer, lr, patience):
    dataset_directory = input("\nINSERIRE IL NOME DELLA CARTELLA IN CUI SI TROVA IL DATASET (sottocartella di "
                              "'Classification/DATASET')").upper()
    labels = int(input("INSERIRE NUMERO DI LABELS PRESENTI NEL DATASET: "))
    cnn_instance = CNN(dataset_directory, model_mk, batch_size, epochs, fl_filter, labels, n_dropout, drop_value,
                       n_layer, lr, patience)
    print("RIEPILOGO RETE:")
    model = cnn_instance.create_model()
    print("VUOI SALVARE I GRAFICI DEI PARAMETRI PER TRAIN/VALIDATION E TEST")
    graph_choice = input("S: SI\nN: NO\n la tua scelta: ").upper()
    print("----------GENERAZIONE BATCH DI DATI ED INIZIO TRAINING----------")
    match graph_choice:
        case "S":
            datagen_list = cnn_instance.datagen()
            history = cnn_instance.train(model, datagen_list[0], datagen_list[1])
            result_plotter = ResultPlotter(history, dataset_directory, cnn_instance.model_filename)
            result_plotter.drawall_graphs()
            score = cnn_instance.test_evaluate(model, datagen_list[2])
            test_plotter = ResultPlotter(score, dataset_directory, cnn_instance.model_filename)
            test_plotter.drawall_test_graphs()

        case "N":
            datagen_list = cnn_instance.datagen()
            cnn_instance.train(model, datagen_list[0], datagen_list[1])


def save_train_val_to_txt(history, epochs, batch_size, lr, k, dataset):
    actual_train_epochs = len(history.history['loss'])
    train_loss = min(history.history['loss'])
    train_accuracy = max(history.history['accuracy'])
    train_precision = max(history.history['precision'])
    train_recall = max(history.history['recall'])
    train_auc = max(history.history['auc'])
    train_f1 = 2 * (train_precision * train_recall) / (train_precision + train_recall)

    val_loss = min(history.history['val_loss'])
    val_accuracy = max(history.history['val_accuracy'])
    val_precision = max(history.history['val_precision'])
    val_recall = max(history.history['val_recall'])
    val_auc = max(history.history['val_auc'])
    val_f1 = 2 * (val_precision * val_recall) / (val_precision + val_recall)

    train_results_path = f"risultati_training_{dataset}_{k}K(epochs={epochs},batch_size={batch_size})"
    #Path(train_results_path).mkdir(parents=True, exist_ok=True)

    training_output_filename = f"{train_results_path}_{lr}_CnnRna.txt"
    with open(training_output_filename, "w") as out_file:
        out_file.write("Risultati training: \n")

    with open(training_output_filename, "a") as out_file:
        out_file.write(
            "epoche: %d\tloss: %.5f\taccuracy: %.5f\tprecision: %.5f\trecall: %.5f\tauroc: %.5f\tF1: %.5f\n" %
            (actual_train_epochs, train_loss, train_accuracy,
             train_precision, train_recall, train_auc, train_f1), )

    with open(training_output_filename, "a") as out_file:
        out_file.write(
            "Risultati Validation:\nepoche:%d\tloss: %.5f\taccuracy: %.5f\tprecision: %.5f\trecall: %.5f\tauroc: %.5f\tF1: %.5f\n" %
            (actual_train_epochs, val_loss, val_accuracy,
             val_precision, val_recall, val_auc, val_f1), )

def save_test_to_txt(score, epochs, batch_size, lr, k, dataset):
    test_loss = score[0]
    test_accuracy = score[1]
    test_precision = score[2]
    test_recall = score[3]
    test_auc = score[4]
    test_f1 = 2 * (test_precision * test_recall) / (test_precision + test_recall)

    test_results_path = f"risultati_test_{dataset}_{k}K(epochs={epochs},batch_size={batch_size})"
    #Path(test_results_path).mkdir(parents=True, exist_ok=True)

    test_output_filename = f"{test_results_path}_{lr}_CnnRna.txt"
    with open(test_output_filename, "w") as out_file:
        out_file.write("Risultati test: \n")

    with open(test_output_filename, "a") as out_file:
        out_file.write("%s\tloss:%.5f\taccuracy:%.5f\tprecision:%.5f\trecall:%.5f\tauroc:%.5f\tF1:%.5f\n" %
                       ("Risultati:", test_loss, test_accuracy,
                        test_precision, test_recall, test_auc, test_f1), )

main()

from scripts.cgr.CGRHandler import CGRHandler
from scripts.cnn.CNN import CNN
from scripts.cnn.ResultPlotter import ResultPlotter

'''
Metodo per l'implementazione della CLI, sfrutta altri metodi per semplificare la lettura
'''


def main():
    model_mk = 1
    batch_size = 32
    epochs = 25
    fl_filter = 32
    n_dropout = 1
    drop_value = 0.5
    n_layer = 3
    lr = 0.0001
    patience = 5

    print("SCEGLI IL PROCESSO DA ESEGUIRE:\n")
    print("1: GENERAZIONE IMMAGINI CGR DA SEQUENZE RNA")
    print("2: CLASSIFICAZIONE IMMAGINI TRAMITE CNN")
    print("3: GENERAZIONE IMMAGINI E CLASSIFICAZIONE")
    process_input = input("la tua scelta: ")
    match process_input:
        case "1":
            imgen_case()
        case "2":
            cnn_case(model_mk, batch_size, epochs, fl_filter, n_dropout, drop_value, n_layer, lr, patience)
        case "3":
            imgen_case()
            cnn_case(model_mk, batch_size, epochs, fl_filter, n_dropout, drop_value, n_layer, lr, patience)


'''
Metodo di supporto a main, contiene il codice per la scelta della generazione delle immagini
'''


def imgen_case():
    fasta_directory = input(
        "\nINSERIRE IL NOME DELLA CARTELLA IN CUI SI TROVANO I FILE FASTA (sottocartella di "
        "'FASTA')").lower()
    images_directory = input(
        "\nINSERIRE IL NOME DELLA CARTELLA IN CUI SI SALVARE LE IMMAGINI (sottocartella di "
        "'IMMAGINI CGR')").lower()
    handler_istance = CGRHandler("RNA", False, False, fasta_directory, images_directory)
    handler_istance.read_files()


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
            history = cnn_instance.train(model, datagen_list[0], datagen_list[1])


main()

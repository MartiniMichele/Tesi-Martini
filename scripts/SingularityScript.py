from scripts.cnn.CNN import CNN
import os
from pathlib import Path
import sys

model_mk = 1
fl_filter = 32
n_dropout = 1
drop_value = 0.5
n_layer = 3
lr = 0.0001
patience = 5


def main(argv):


    dataset_directory = "C:/Users/Michele/Documents/GitHub/Tesi-Martini/Classification/LSU_DATASET"
    labels = 3
    epochs = 1
    batch_size = 2048

    cnn_instance = CNN(dataset_directory, model_mk, batch_size, epochs, fl_filter, labels, n_dropout, drop_value,
                       n_layer, lr, patience)
    model = cnn_instance.create_singularity_model()

    datagen_list = cnn_instance.datagen()
    history = cnn_instance.train(model, datagen_list[0], datagen_list[1])
    score = cnn_instance.test_evaluate(model, datagen_list[2])

    ####################################################################################################
    ####################################################################################################
    #
    #                             salvataggio risultati training/validation
    #
    ####################################################################################################
    ####################################################################################################
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



    train_results_path = f"risultati_training(epochs={epochs},batch_size={batch_size})"
    Path(train_results_path).mkdir(parents=True, exist_ok=True)

    training_output_filename = f"{train_results_path}_{lr}_CnnRna.txt"
    with open(training_output_filename, "w") as out_file:
        out_file.write("Risultati training: \n")

    with open(training_output_filename, "a") as out_file:
        out_file.write("epoche: %d\tloss: %.5f\taccuracy: %.5f\tprecision: %.5f\trecall: %.5f\tauroc: %.5f\tF1: %.5f\n" %
                       (actual_train_epochs, train_loss, train_accuracy,
                        train_precision, train_recall, train_auc, train_f1), )

    with open(training_output_filename, "a") as out_file:
        out_file.write("Risultati Validation:\nepoche:%d\tloss: %.5f\taccuracy: %.5f\tprecision: %.5f\trecall: %.5f\tauroc: %.5f\tF1: %.5f\n" %
                       (actual_train_epochs, val_loss, val_accuracy,
                        val_precision, val_recall, val_auc, val_f1), )

    ####################################################################################################
    ####################################################################################################
    #
    #                                   salvataggio risultati test
    #
    ####################################################################################################
    ####################################################################################################
    test_loss = score[0]
    test_accuracy = score[1]
    test_precision = score[2]
    test_recall = score[3]
    test_auc = score[4]
    test_f1 = 2 * (test_precision * test_recall) / (test_precision + test_recall)


    test_results_path = f"risultati_test(epochs={epochs},batch_size={batch_size})"
    Path(test_results_path).mkdir(parents=True, exist_ok=True)

    test_output_filename = f"{test_results_path}_{lr}_CnnRna.txt"
    with open(test_output_filename, "w") as out_file:
        out_file.write("Risultati test: \n")

    with open(test_output_filename, "a") as out_file:
        out_file.write("%s\tloss:%.5f\taccuracy:%.5f\tprecision:%.5f\trecall:%.5f\tauroc:%.5f\tF1:%.5f\n" %
                       ("Risultati:", test_loss, test_accuracy,
                        test_precision, test_recall, test_auc, test_f1), )


if __name__ == "__main__":
    main(sys.argv[1:])

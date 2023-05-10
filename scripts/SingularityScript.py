from scripts.cnn.CNN import CNN
import os
from pathlib import Path

model_mk = 1
batch_size = 32
epochs = 25
fl_filter = 32
n_dropout = 1
drop_value = 0.5
n_layer = 3
lr = 0.0001
patience = 5

source_path = Path(__file__).resolve()
source_dir = Path(source_path.parent.parent)
dataset_directory = Path(str(source_dir) + "/DATASET")
labels = 3

cnn_instance = CNN(dataset_directory, model_mk, batch_size, epochs, fl_filter, labels, n_dropout, drop_value,
                   n_layer, lr, patience)
model = cnn_instance.create_model()
datagen_list = cnn_instance.datagen()
history = cnn_instance.train(model, datagen_list[0], datagen_list[1])
score = cnn_instance.test_evaluate(model, datagen_list[2])

print(f"risultati training {history}")
print(f"risultati test {score}")
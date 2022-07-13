

import random
from argparse import ArgumentParser
import multiprocessing

import torch
import torchio as tio

import pytorch_lightning as pl
from MPLearn.virtual_staining.data_module import VirtualStainingDataModule
from MPLearn.virtual_staining.linear_model import Linear2p5

seed = 48104
random.seed(seed)
torch.manual_seed(seed)

# dimensions
# dim_x = 2160
# dim_y = 2560
# dpc_dim_z = 30
# stain_dim_z = 10

class Crapifier:
    def __init__(self):
        self.stain_transform = tio.RescaleIntensity(
            out_min_max = (0, 1))
    def __call__(self, dpc_subject, stain_subject):
        return (
            dpc_subject,
            self.stain_transform(stain_subject))


data_module = VirtualStainingDataModule(
    dpc_images_dir = "raw_data/virtual_staining/20200530T160342_48_hour_merged_well4_color1",
    stain_images_dir = "raw_data/virtual_staining/20200530T160342_48_hour_merged_well4_color1",
    queue_length = 256,
    queue_num_workers = 10,
    samples_per_volume = 128,
    patch_dim_x = 25,
    patch_dim_y = 25,            
    batch_size = 16,
    transform = Crapifier())
    
data_module.setup(stage='train')

root_dir = 'runs'
parser = ArgumentParser(add_help=False)
parser = Linear2p5.add_model_specific_args(parser)
parser = pl.Trainer.add_argparse_args(parser)
hparams = parser.parse_args(
    args=[
        "--kernel_size", "3",
        "--learning_rate", "0.002"])

model = Linear2p5(hparams)

trainer = pl.Trainer.from_argparse_args(
    hparams,
    gpus = 1,
    log_every_n_steps=1)

trainer.fit(model, datamodule=data_module)

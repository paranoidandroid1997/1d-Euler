import torch.nn.functional as F
import lightning.pytorch as pl
import torch
from torch import optim, nn

class Shock_Cell_Classifier(pl.LightningModule):
    def __init__(self, input_dim, lr = 1e-4, wd = 1.0e-5):
        super().__init__()
        self.save_hyperparameters()

        # self.fc1 = nn.Linear(self.hparams.input_dim, 16)
        # self.fc2 = nn.Linear(16, 2)
        
        self.fc1 = nn.Linear(self.hparams.input_dim, 8)
        self.fc2 = nn.Linear(8, 2)

        self.dropout5 = nn.Dropout(0.50)
        self.dropout2 = nn.Dropout(0.20)

        nn.init.xavier_uniform_(self.fc1.weight)
        nn.init.xavier_uniform_(self.fc2.weight)
        
        #Initialize bias terms to zero
        nn.init.constant_(self.fc1.bias, 0)
        nn.init.constant_(self.fc2.bias, 0)

        self.loss_fn = nn.CrossEntropyLoss()

        # self.logger.log_hyperparams()
    def forward(self, x):
        # x = self.dropout2(x)
        x = self.fc1(x)
        x = F.leaky_relu(x)
        #x = self.dropout2(x)
        x = self.fc2(x)

        return x

    def training_step(self, batch, batch_idx):
        # training_step defines the train loop.
        # it is independent of forward
        inputs, labels = batch
        
        output = self.forward(inputs)

        
        loss = self.loss_fn(output, labels.long())
        # # Logging to TensorBoard (if installed) by default
        self.log("train_loss", loss)
        return loss

    def validation_step(self, batch, batch_idx):
        # this is the validation loop
        inputs, labels = batch
        
        output = self.forward(inputs)
        
        val_loss = self.loss_fn(output, labels.long())
        self.log("val_loss", val_loss)

    def configure_optimizers(self):
        optimizer = optim.Adam(self.parameters(), lr=self.hparams.lr, weight_decay=self.hparams.wd)
        return optimizer
    
    def initialize():
        stencil = 5
        classifier = Shock_Cell_Classifier(stencil * 4)
        classifier.load_state_dict(torch.load("./nn-02.pt"))
        classifier.eval()

        return classifier
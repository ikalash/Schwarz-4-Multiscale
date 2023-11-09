import argparse

from pinnschwarz.trainer import Trainer


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="PINN-Schwarz", description="Schwarz-based training of PINN-PINN coupling")

    parser.add_argument("parameter_file")
    parser.add_argument("outdir")
    args = parser.parse_args()

    train_mod = Trainer(args.parameter_file, args.outdir, "hyper.yaml", False)
    train_mod.train()

"""Console script for tf_circuit_driver."""
import click
import sys
from tf_circuit_driver.scripts.bootstrap_db import bootstrap, build_protein_nodes
from tf_circuit_driver.scripts.calculate import calculate


@click.command()
@click.option('--actions', help="filepath to the StringDB proteins actions table")
@click.option('--mappings', help="filepath to the StringDB protein mapping")
@click.option('--targets', help="directory path to your unzipped targets directory form the MARA report")
@click.option('--dge', help="filepath to your differential expression file")
@click.option('--output', help="filepath where your final results will be")
@click.option('--mode', help="Running MODE", type=click.Choice(['init', 'calculate']))
def main(actions, mappings, targets, mode, dge, output):
    if mode == 'init':
        bootstrap(actions, mappings, targets, dge)
    elif mode == 'calculate':
        calculate(output)


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover

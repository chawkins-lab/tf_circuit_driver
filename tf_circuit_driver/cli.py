"""Console script for tf_circuit_driver."""
import click
import sys
from tf_circuit_driver.scripts.bootstrap_db import bootstrap_ppi, parse_dge, construct_mara
from tf_circuit_driver.scripts.calculate import calculate


@click.group()
def main():
    pass


@click.command()
@click.option('--actions', help="filepath to the StringDB proteins actions table")
@click.option('--mappings', help="filepath to the StringDB protein mapping")
def init_ppi(actions, mappings):
    click.echo('Initializing Database')
    bootstrap_ppi(actions, mappings)


@click.command()
@click.option('--targets', help="directory path to your unzipped targets directory form the MARA report")
def init_mara(targets):
    click.echo('Mapping MARA targets directory')
    construct_mara(targets)


@click.command()
@click.option('--dge', help="filepath to your differential expression file")
def map_dge(dge):
    click.echo('Mapping Differential Expression Scores')
    parse_dge(dge_filepath=dge)


@click.command()
@click.option('--output', help="filepath to write output results")
def calculate_scores(output):
    calculate(output)


main.add_command(init_ppi)
main.add_command(init_mara)
main.add_command(map_dge)
main.add_command(calculate_scores)
if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover

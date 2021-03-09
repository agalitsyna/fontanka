import click
from . import cli

@cli.command()
@click.argument(
    "cool_path",
    metavar="COOL_PATH",
    type=str
)
@click.option(
    "--kernel-size",
    help="Kernel size",
    type=int,
)
def call_fountains(
        cool_path,
        kernel_size
        ):
    print(cool_path)

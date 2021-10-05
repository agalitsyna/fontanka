import click
from .. import __version__
from .._logging import get_logger

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}


@click.version_option(__version__, "-V", "--version")
@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """Type -h or --help after subcommand."""
    pass


from . import call_fountains

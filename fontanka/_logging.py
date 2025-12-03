# Based on ipython traitlets
import logging
import sys

logging.basicConfig(
    stream=sys.stdout, 
    level=logging.INFO,
    format="[%(asctime)s] %(name)s:%(levelname)s: %(message)s",
    datefmt='%H:%M:%S'
)

_loggers = {}


def get_logger(name="fontanka"):
    global _loggers
    if name not in _loggers:
        _loggers[name] = logging.getLogger(name)
        # Add a NullHandler to silence warnings about not being
        # initialized.
        _loggers[name].addHandler(logging.NullHandler())
    return _loggers[name]

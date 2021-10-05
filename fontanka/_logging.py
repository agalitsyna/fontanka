# Based on ipython traitlets
import logging

_loggers = {}


def get_logger(name="fontanka"):
    global _loggers
    if name not in _loggers:
        _loggers[name] = logging.getLogger(name)
        # Add a NullHandler to silence warnings about not being
        # initialized, per best practice for libraries.
        _loggers[name].addHandler(logging.NullHandler())
    return _loggers[name]

from setuptools import setup
import re
import os
import io


def _read(*parts, **kwargs):
    filepath = os.path.join(os.path.dirname(__file__), *parts)
    encoding = kwargs.pop("encoding", "utf-8")
    with io.open(filepath, encoding=encoding) as fh:
        text = fh.read()
    return text


def get_version():
    version = re.search(
        r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
        _read("fontanka", "__init__.py"),
        re.MULTILINE,
    ).group(1)
    return version


setup(
    name="fontanka",
    version=get_version(),
    author="agalicina",
    license="MIT",
    py_modules=["fontanka"],
    install_requires=[
        "Click",
    ],
    entry_points={
        "console_scripts": [
            "fontanka = fontanka.cli:cli",
        ]
    },
)

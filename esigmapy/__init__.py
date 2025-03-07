from __future__ import absolute_import

from . import blend, generator, legacy, utils
from .generator import *


def get_version_information():
    import os

    version_file = os.path.join(
        os.path.dirname(os.path.dirname(__file__)), "esigmapy/.version"
    )
    try:
        with open(version_file, "r") as f:
            return f.readline().rstrip()
    except EnvironmentError:
        print("No version information file '.version' found")


__version__ = get_version_information()

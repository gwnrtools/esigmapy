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


# pycbc wrapper
import inspect
_params = inspect.signature(get_imr_esigma_waveform).parameters.keys()
def pycbc_esigma(**params):
    pt = {p: params[p] for p in _params if p in params}
    pt['f_ref'] = None

    gen_wav = get_imr_esigma_waveform
    if 'skip_merger' in params:
        gen_wav = get_inspiral_esigma_waveform

    return gen_wav(pt.pop('mass1'), pt.pop('mass2'),
                            pt.pop('f_lower'), pt.pop('delta_t'),
                            **pt, modes_to_use=[(2,2), (2,-2)])

__version__ = get_version_information()

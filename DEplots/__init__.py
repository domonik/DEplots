from .diffexp import *
import sys
from . import diffexp as _newname
sys.modules[__name__ + '.volcano'] = _newname
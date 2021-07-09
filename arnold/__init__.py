# check if SAENO path is stored
import os
from distutils.sysconfig import get_python_lib

cfg_file = get_python_lib() + '/saeno.cfg'

if os.path.exists(cfg_file):
    with open(cfg_file, 'r') as f:
        SAENOPATH = f.readlines()[0].strip()
else:
    SAENOPATH = '.'

# import package files
from . import experiment
from . import force
from . import materials
from . import mesh
from . import simulation
from . import simulation_saenopy
from . import utils
from . import image
from . import regularization
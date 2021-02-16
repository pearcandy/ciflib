__copyright__ = 'Copyright(c) 2020 PearCandy'
__version__ = '0.1.2'
__lisence__ = 'MIT'
__author__ = 'PearCandy'
__author_email__ = 'y.nishi1980@gmail.com'
__all__ = ['src']

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

from .ciflib import main

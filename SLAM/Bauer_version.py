# My best shot at getting SLAM! to run 

# it earned the explaimation mark after the third album, duh

import numpy as np 
import matplotlib.pyplot as plt 

from the_model import *
from make_forcing import * 

# it seems like write_netCDF makes a CDF of noise, that we can
# feed to SLAM to process 

# here it goes

years = 10 

write_netCDF(years)


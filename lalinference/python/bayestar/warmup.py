import threading
import numpy as np
from lalinference.bayestar.sky_map import toa_phoa_snr

thread = threading.Thread(
    target=toa_phoa_snr,
    args=(0., [], [], [], [], [], [],
    np.empty((0, 3, 3), dtype=np.float32),
    np.empty((0, 3)), [], 0., 1., 2, 1))
thread.start()

finish = thread.join


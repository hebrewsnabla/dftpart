from dftpart import scfeda
import numpy as np
from QCKit import logger
                     
subeda = scfeda.EDA()
subeda.method = ["hf", "6-31gs"]
subeda.gjf = '32ane_subsys/32ane_5.gjf'
subeda.output =  subeda.gjf[:-4]
subatm_E, subE, conv = subeda.kernel()
atom_E = np.zeros(98)
resultlog = '32ane-result.log'
with open(resultlog,'a') as g:
    logger.slog(g, "center atom energies:")
    cen_intot = [21, 22, 23, 24, 76, 77, 78, 79, 80, 81, 82, 83]
    cen_insub = [5, 6, 7, 8, 24, 25, 26, 27, 28, 29, 30, 31]
    for i in range(len(cen_intot)):                        
        atom_E[cen_intot[i]-1] = subatm_E[cen_insub[i]-1]                        
        logger.slog(g, "%d %.10f", cen_intot[i], atom_E[cen_intot[i]-1])

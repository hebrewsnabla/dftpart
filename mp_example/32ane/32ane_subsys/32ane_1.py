from dftpart import scfeda
import numpy as np
from QCKit import logger
                     
subeda = scfeda.EDA()
subeda.method = ["hf", "6-311gss"]
subeda.gjf = '32ane_subsys/32ane_1.gjf'
subeda.output =  subeda.gjf[:-4]
subatm_E, subE, conv = subeda.kernel()
atom_E = np.zeros(98)
resultlog = '32ane-result_1.log'
with open(resultlog,'a') as g:
    logger.slog(g, "## center atom energies in subsystem 1:")
    cen_intot = [1, 2, 3, 4, 5, 6, 7, 8, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53]
    cen_insub = [1, 2, 3, 4, 5, 6, 7, 8, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
    for i in range(len(cen_intot)):                        
        atom_E[cen_intot[i]-1] = subatm_E[cen_insub[i]-1]                        
        logger.slog(g, "%d %.10f", cen_intot[i], atom_E[cen_intot[i]-1])

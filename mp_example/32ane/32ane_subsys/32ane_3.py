from dftpart import scfeda
import numpy as np
from QCKit import logger
                     
subeda = scfeda.EDA()
subeda.method = ["hf", "6-311gss"]
subeda.gjf = '32ane_subsys/32ane_3.gjf'
subeda.output =  subeda.gjf[:-4]
subatm_E, subE, conv = subeda.kernel()
atom_E = np.zeros(98)
resultlog = '32ane-result_3.log'
with open(resultlog,'a') as g:
    logger.slog(g, "## center atom energies in subsystem 3:")
    cen_intot = [13, 14, 15, 16, 36, 62, 63, 64, 65, 66, 67, 68]
    cen_insub = [5, 6, 7, 8, 16, 25, 26, 27, 28, 29, 30, 31]
    for i in range(len(cen_intot)):                        
        atom_E[cen_intot[i]-1] = subatm_E[cen_insub[i]-1]                        
        logger.slog(g, "%d %.10f", cen_intot[i], atom_E[cen_intot[i]-1])

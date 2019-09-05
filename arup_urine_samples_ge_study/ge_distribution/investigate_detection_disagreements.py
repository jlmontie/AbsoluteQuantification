import json
import numpy as np
# Used to build detection_disagreements.txt
sample_summary = '/uufs/chpc.utah.edu/common/home/u0002613/AbsoluteQuantification/AbsoluteQuantification/data/arup_urine_summary_files_with_quantification/190828-1-1-IDBD-D100389-d-13-AHMW75AFXY-CGTTATTCTA-AACCTTATGG-5d6830eb-r.rna.bacterial.dxsm.out.summary'
idbd_detection_taxid = 562
with open(sample_summary) as file:
    for line in file:
        obj = json.loads(line)
        if obj['taxid'] == idbd_detection_taxid:
            print(obj['absolute_quant'])
            print(np.log10(obj['absolute_quant']))
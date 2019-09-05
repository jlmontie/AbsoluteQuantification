import json
import numpy as np
# Used to build detection_disagreements.txt
sample_summary = '/Users/jmontgomery/Desktop/tmp_summary_quant/190903-1-1-IDBD-D100414-d-04-AHMV3HAFXY-TCATAGATTG-CACCTTAATC-5d700753-r.rna.bacterial.dxsm.out.summary'
idbd_detection_taxid = 287
with open(sample_summary) as file:
    for line in file:
        obj = json.loads(line)
        # print(obj.keys())
        # print(obj['median_depth'])
        # print(obj['coverage'])
        # break
        if obj['taxid'] == idbd_detection_taxid:
            print(obj['absolute_quant'])
            print(np.log10(obj['absolute_quant']))
import json

bac_summary = '/uufs/chpc.utah.edu/common/home/u0002613/synergy_validation/190828-1-1-IDBD-D100376-d-01-AHMW75AFXY-TGATTATACG-TGTAATCGAC-5d6e759f.rna.bacterial.dxsm.out.summary'
with open(bac_summary) as file:
    for line in file:
        obj = json.loads(line)
        if obj['taxid'] == 562:
            print(f"name: {obj['name']}\nabsolute_quant: {obj['absolute_quant']}\nabsolute_quant_bin: {obj['absolute_quant_bin']}")
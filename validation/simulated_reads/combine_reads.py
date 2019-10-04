import numpy as np

files = np.array([
    'T7_090918RESP115-18250402467-d_NC_001604.fa',  # 0 - IC
    '562_e_coli_quant_6.5_depth_4.fa',  # 1 - Bacteria A
    '562_e_coli_quant_7.5_depth_34.fa',  # 2 - Bacteria A
    '562_e_coli_quant_8.5_depth_329.fa',  # 3 - Bacteria A
    '573_k_pneumoniae_quant_6.5_depth_4.fa',  # 4 - Bacteria B
    '573_k_pneumoniae_quant_7.5_depth_39.fa',  # 5 - Bacteria B
    '573_k_pneumoniae_quant_8.5_depth_376.fa',  # 6 - Bacteria B
    '5476_c_albicans_quant_6.5_depth_24.fa',  # 7 - Fungi C
    '5476_c_albicans_quant_7.5_depth_235.fa',  # 8 - Fungi C
    '5476_c_albicans_quant_8.5_depth_2278.fa',  # 9 - Fungi C
    '5482_c_tropicalis_quant_6.5_depth_24.fa',  # 10 - Fungi D
    '5482_c_tropicalis_quant_7.5_depth_235.fa',  # 11 - Fungi D
    '5482_c_tropicalis_quant_8.5_depth_2278.fa'  # 12 - Fungi D
])

combo_idx = [
    [0],  # Q-IC
    [0, 2],  # Q-Bacteria-1
    [0, 5],  # Q-Bacteria-2
    [0, 2, 5],  # Q-Bacteria-Mid
    [0, 3, 4],  # Q-Bacteria-HiLo-1
    [0, 1, 6],  # Q-Bacteria-HiLo-2
    [0, 8],  # Q-Fungi-1
    [0, 11],  # Q-Fungi-2
    [0, 8, 11],  # Q-Fungi-Mid
    [0, 9, 10],  # Q-Fungi-HiLo-1
    [0, 7, 12],  # Q-Fungi-HiLo-2
    [0, 2, 5, 8, 11],  # Q-Mixed-Mid
    [0, 3, 4, 9, 10],  # Q-Mixed-HiLo-1
    [0, 1, 6, 7, 12],  # Q-Mixed-HiLo-2
    [2, 5, 8, 11]  # Q-Mixed-Mid-No-IC
]

combo_names = np.array([
    'Q-IC.fa',
    'Q-Bacteria-1.fa',
    'Q-Bacteria-2.fa',
    'Q-Bacteria-Mid.fa',
    'Q-Bacteria-HiLo-1.fa',
    'Q-Bacteria-HiLo-2.fa',
    'Q-Fungi-1.fa',
    'Q-Fungi-2.fa',
    'Q-Fungi-Mid.fa',
    'Q-Fungi-HiLo-1.fa',
    'Q-Fungi-HiLo-2.fa',
    'Q-Mixed-Mid.fa',
    'Q-Mixed-HiLo-1.fa',
    'Q-Mixed-HiLo-2.fa',
    'Q-Mixed-Mid-No-IC.fa'
])

for idx_ls, name in zip(combo_idx, combo_names):
    with open(name, 'w') as file_out:
        for idx in idx_ls:
            with open(files[idx]) as file_in:
                for line in file_in:
                    file_out.write(line)

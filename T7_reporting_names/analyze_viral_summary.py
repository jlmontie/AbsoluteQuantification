import json
total_read_count = 1e6
with open('T7_090918RESP115-18250402467-d_NC_001604.rna.viral.dxsm.out.summary') as file, \
    open('t7_simulate_reads_summary_info.txt', 'w') as outfile:
    outfile.write('taxid\tncbi_name\treporting_id\tcompound_id\tread_count\tread_pcnt\n')
    for line in file:
        org_info = json.loads(line)
        write_str = f"{org_info['taxid']}\t{org_info['ncbi_name']}\t{org_info['reporting_id']}\t{org_info['compound_id']}\t{org_info['read_count']}\t{100 * org_info['read_count'] / total_read_count:0.1f}\n"
        outfile.write(write_str)


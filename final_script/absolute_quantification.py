import os
print(taxid_reporting_name_info['26706_10760'])


def absolute_quantification(config, summary_paths, ic_reporting_ids):
    # load taxid reporting name info
    taxid_reporting_name_info_path = os.path.join(
        config['top_level_dir'], config['taxid_reporting_names'])
    # {"compound_id"}->{taxid:, reporting_id:, reporting_name:, class_type:, subclass:}
    taxid_reporting_name_info = {}
    with open(taxid_reporting_name_info_path) as f:
        f.readline()  # skip header
        for line in f:
            data = line.strip().split("\t")
            for compound_id in data[3].split(","):
                taxid_reporting_name_info[compound_id] = {
                    "taxid": data[2],
                    "reporting_id": data[1],
                    "reporting_name": data[0],
                    "class_type": data[4],
                    "subclass": [x for x in data[5].split(",")],
                    "nucleic_acid": data[6],
                    "RNA": data[7],
                    "DNA": data[8]
                }

                # overload viral dbs tag
                if taxid_reporting_name_info[compound_id]["class_type"] == "viral":
                    taxid_reporting_name_info[compound_id]["DNA"] = taxid_reporting_name_info[compound_id]["RNA"]

                if taxid_reporting_name_info[compound_id]["RNA"] == ".":
                    del taxid_reporting_name_info[compound_id]["RNA"]
                if taxid_reporting_name_info[compound_id]["DNA"] == ".":
                    del taxid_reporting_name_info[compound_id]["DNA"]
    f.close()

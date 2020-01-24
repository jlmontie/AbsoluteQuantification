# banks server
SCRIPT=/opt/classification_synergy/utils/explify_v2/exbox_pipeline_v2.py
CONFIG=/home/jmontgomery/AbsoluteQuantification/validation/pipeline/synergy_validation_jan_2020/config_banks.json
YAML=/home/jmontgomery/AbsoluteQuantification/validation/pipeline/synergy_validation_jan_2020/input_in_silico_validation.yaml

python $SCRIPT $YAML $CONFIG
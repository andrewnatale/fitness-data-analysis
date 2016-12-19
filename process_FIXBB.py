import sys
import os
from itertools import chain
from StateDataTools import *
import heatmaps

data_path = '/Users/anatale/Documents/school/UCSF/Kortemme_lab/data/second_run_scores'

yeast_gsp1_seq = '\
MSAPAANGEVPTFKLVLVGDGGTGKTTFVKRHLTGEFEKKYIATIGVEVHPLSFYTNFGE\
IKFDVWDTAGQEKFGGLRDGYYINAQCAIIMFDVTSRITYKNVPNWHRDLVRVCENIPIV\
LCGNKVDVKERKVKAKTITFHRKKNLQYYDISAKSNYNFEKPFLWLARKLAGNPQLEFVA\
SPALAPPEVQVDEQLMQQYQQEMEQATALPLPDEDDADL'

position_ids = [101,102,103,105,106,107,108,109,122,123,124,125,126,127]

state_names = ['apo_GEF','GTP_GAP','GDP','GTP_importin']
#state_names = ['apo_GEF',]
processor = StateProcessor(position_ids)
processor.seq_from_string(yeast_gsp1_seq, 'yeast GSP1')
processor.generate_labels()

for idx,state in enumerate(state_names):
    processor.load_state_pickle(os.path.join(data_path, 'state_%d_data_v3.pkl' % (idx+1)), state_name=state)

#processor.write_microstate_seqs()

#processor.write_macrostate_energy_tsv()
#processor.write_macrostate_energy_tsv(outname='testing.tsv', position_range=range(122,128,1))
processor.write_microstate_energy_tsv(outname='testing_microstates.tsv', position_range=range(122,128,1))

# for state in processor.state_list:
#     print(state.state_name)
#     heatmaps.single_map(
#     state.avg_state_array,
#     processor.mutation_labels,
#     processor.sequence_labels,
#     use_mpoint=True,
#     set_limits=(20,-5),
#     show_cbar=True
#     )
#
# for state in processor.state_list:
#     print(state.state_name)
#     heatmaps.single_map(
#     state.exp_state_array,
#     processor.mutation_labels,
#     processor.sequence_labels,
#     use_mpoint=True,
#     mpoint=1.0,
#     set_limits=(2,0),
#     show_cbar=True
#     )

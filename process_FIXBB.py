import sys
import os
from itertools import chain
from StateDataTools import *
import heatmaps

output_dir = '/Users/anatale/school/UCSF/Kortemme_lab/code/multi-state-design/opt_rnd3'

data_path = '/Users/anatale/school/UCSF/Kortemme_lab/data/combo_scores/new_final'

yeast_gsp1_seq = '\
MSAPAANGEVPTFKLVLVGDGGTGKTTFVKRHLTGEFEKKYIATIGVEVHPLSFYTNFGE\
IKFDVWDTAGQEKFGGLRDGYYINAQCAIIMFDVTSRITYKNVPNWHRDLVRVCENIPIV\
LCGNKVDVKERKVKAKTITFHRKKNLQYYDISAKSNYNFEKPFLWLARKLAGNPQLEFVA\
SPALAPPEVQVDEQLMQQYQQEMEQATALPLPDEDDADL'

position_ids = [101,102,103,104,105,106,107,108,109,122,123,124,125,126,127]

state_names = ['1i2m', '1a2k', '1k5d', '3gj0', 'importin', 'composite']
# initialize data processor
processor = StateProcessor(position_ids)
processor.seq_from_string(yeast_gsp1_seq, 'yeast GSP1')
processor.generate_labels()

# # load data from all states
# for idx,state in enumerate(state_names):
#     processor.load_state_pickle(os.path.join(data_path, 'state_%d_data_v5.pkl' % (idx+1)), state_name=state)

processor.load_state_pickle(os.path.join(data_path, 'state_3_data_v5.pkl'), state_name='Ran-GTP-GAP')
processor.load_state_pickle(os.path.join(data_path, 'state_4_data_v5.pkl'), state_name='Ran-GDP')
processor.load_state_pickle(os.path.join(data_path, 'state_2_data_v5.pkl'), state_name='Ran-GDP-NTF2')
processor.load_state_pickle(os.path.join(data_path, 'state_1_data_v5.pkl'), state_name='Ran-GEF')
processor.load_state_pickle(os.path.join(data_path, 'state_5_data_v5.pkl'), state_name='Ran-GTP-karyopherin')
processor.load_state_pickle(os.path.join(data_path, 'state_6_data_v5.pkl'), state_name='Composite')

# for state in processor.state_list:
#     print(state.state_name, state.microstate_count)
#     for microstate in state.microstate_list:
#         print(microstate.substate_label, microstate.us_name, microstate.find_best_seq())

def gen_all_files():
    current_dir = os.getcwd()
    os.chdir(output_dir)

    processor.write_microstate_seqs(position_range=range(122,128,1), label_string='122_127')
    processor.write_microstate_seqs(position_range=range(101,110,1), label_string='101_109')
    processor.write_microstate_seqs(position_range='all', label_string='all')

    processor.write_macrostate_energy_tsv(outname='macro_energy_122_127.tsv', position_range=range(122,128,1))
    processor.write_macrostate_energy_tsv(outname='macro_energy_101_109.tsv', position_range=range(101,110,1))
    processor.write_macrostate_energy_tsv(outname='macro_energy_all.tsv', position_range='all')

    processor.write_microstate_energy_tsv(outname='micro_energy_122_127.tsv', position_range=range(122,128,1))
    processor.write_microstate_energy_tsv(outname='micro_energy_101_109.tsv', position_range=range(101,110,1))
    processor.write_microstate_energy_tsv(outname='micro_energy_all.tsv', position_range='all')

    os.chdir(current_dir)

def plot_heatmaps():
    # avg
    heatmaps.new_multi_map(
    [state.avg_state_array[:,:9] for state in processor.state_list],
    processor.mutation_labels,
    processor.sequence_labels[:9],
    ['a','b','c','d','e','f'],
    use_mpoint=True,
    mpoint=0,
    set_limits=(20,-5),
    show_cbar=True
    )
    # exp scale
    heatmaps.new_multi_map(
    [state.exp_state_array[:,9:] for state in processor.state_list],
    processor.mutation_labels,
    processor.sequence_labels[9:],
    [state.state_name for state in processor.state_list],
    use_mpoint=True,
    mpoint=1.0,
    set_limits=(2,0),
    show_cbar=True,
    reverse_cmap=True
    )

#gen_all_files()
plot_heatmaps()
#processor.write_microstate_energy_tsv(outname='micro_energy_all_remap.tsv', position_range='all')

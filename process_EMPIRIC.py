from EMPIRICtools import *
import heatmaps

native_seq_file = '../../notebook/logo/P32835.fasta'
fitness_data_file = '../../notebook/logo/gsp_log_data_101-140.csv'

processor = EMPIRICProcessor()
processor.seq_from_fasta(native_seq_file)
processor.load_data_csv(fitness_data_file)
processor.generate_labels()

processor.rescale()
processor.generate_freqencies()

processor.errorstats()
# print(processor.wt_vals.shape)
# print(processor.stop_vals.shape)
# print('wt mean', processor.wt_mean)
# print('wt std', processor.wt_std)
# print('stop mean', processor.stop_mean)
# print('stop std', processor.stop_std)
#processor.pssm_out(position_range=(122,127), outfile='122-127_pssm_out.txt')
# processor.pssm_out()
#processor.fasta_out(scale_factor=1000, position_range=(122,127), outfile='EMPIRIC_freq_122_127.fasta')
# processor.fasta_out(scale_factor=1000, position_range=(101,109), outfile='EMPIRIC_freq_101_109.fasta')
#processor.fasta_out(scale_factor=1000, position_range=[101,102,103,104,105,106,107,108,109,122,123,124,125,126,127], outfile='EMPIRIC_freq_101_127.fasta')
#processor.fasta_out(scale_factor=80, position_range=None, outfile='EMPIRIC_all_80.fasta')

# heatmaps.single_map(processor.rescaled_array[1:,21:27],
#                    processor.mutation_labels[1:],
#                    processor.sequence_labels[21:27],
#                    use_mpoint=True,
#                    mpoint=1.0,
#                    set_limits=None,
#                    show_cbar=True,
#                    reverse_cmap=True)

# heatmaps.single_map(processor.rescaled_array[1:,0:9],
#                    processor.mutation_labels[1:],
#                    processor.sequence_labels[0:9],
#                    use_mpoint=True,
#                    mpoint=1.0,
#                    set_limits=None,
#                    show_cbar=True,
#                    reverse_cmap=True)

#
# mydata.rescale()
#
# heatmap.single_map(processor.rescaled_array[1:,:],
#                    processor.mutation_labels[:-1],
#                    processor.sequence_labels,
#                    show_cbar=True)
#
# mydata.generate_freqencies()
#
# heatmap.single_map(processor.frequency_array[1:,:],
#                    processor.mutation_labels[:-1],
#                    processor.sequence_labels,
#                    show_cbar=True)

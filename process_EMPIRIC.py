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

#processor.pssm_out(position_range=(122,127), outfile='122-127_pssm_out.txt')
# processor.pssm_out()
processor.fasta_out(scale_factor=1000, position_range=(122,127), outfile='EMPIRIC_freq_122_127.fasta')

# heatmaps.single_map(processor.fitness_array,
#                    processor.mutation_labels,
#                    processor.sequence_labels,
#                    use_mpoint=True,
#                    show_cbar=True)
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

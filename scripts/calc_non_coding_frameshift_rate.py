import argparse
import csv
import pandas as pd


def get_frameshift_info(fs_df, bins):
    fs_df['indel len'] = fs_df['End_Position'] - fs_df['Start_Position']

    # count the number INDELs with length non-dividable by 3
    num_indels = []
    indel_len = []
    num_categories = 0
    i = 1
    while(num_categories<bins):
        # not inframe length indel
        if i%3:
            if num_categories != bins-1:
                tmp_num = len(fs_df[fs_df['indel len']==i])
            else:
                tmp_num = len(fs_df[(fs_df['indel len']>=i) & ((fs_df['indel len']%3)>0)])
            num_indels.append(tmp_num)
            indel_len.append(i)
            num_categories += 1
        i += 1
    return indel_len, num_indels


def get_genome_length(path):
    genome_len = 0
    with open(path) as handle:
        for line in csv.reader(handle, delimiter='\t'):
            if len(line) > 1:
                chrom_len = int(line[1])
                genome_len += chrom_len
    return genome_len


def get_black_list_length(path):
    black_list_len = 0
    with open(path) as handle:
        for line in csv.reader(handle, delimiter='\t'):
            tmp_len = int(line[2]) - int(line[1])
            black_list_len += tmp_len
    return black_list_len


def parse_arguments():
    info = 'Calculates background non-coding frameshift rate stratified by length.'
    parser = argparse.ArgumentParser(description=info)
    help_str = ('BED file containing non-overlapping regions where INDELs are '
                'not counted as non-coding (i.e. gene annotations, etc.)')
    parser.add_argument('-b', '--black-list',
                        type=str, required=True,
                        help=help_str)
    help_str = 'File containing chromosome lengths in base pairs'
    parser.add_argument('-g', '--genome-file',
                        type=str, required=True,
                        help=help_str)
    help_str = ('Expected average coverage for genome sequencing. This should'
               'be sufficient coverage to make calls for indels. Thus the '
               'number of bases with >30x may be reasonable with adjustments '
                'for sample purity. (Default: 1.0)')
    parser.add_argument('-c', '--coverage',
                        type=float, default=1.0,
                        help=help_str)
    help_str = 'Tab-delimited Mutations file containing non-coding indels'
    parser.add_argument('-i', '--indels',
                        type=str, required=True,
                        help=help_str)
    help_str = ('Minimum threshold of "non-coding" indels for a sample to be used '
                'in background rate calculations. This prevents issues with gene '
                'annotations spurious classifying an entire sample as containing '
                'non coding indels.')
    parser.add_argument('-t', '--threshold',
                        type=int, default=10,
                        help=help_str)
    help_str = ('Number of bins to stratify frameshift length by. Each bin '
                'encompasses a single frameshift length (e.g. length=4) up until '
                'the last bin which encorporates every length equivalent and above.')
    parser.add_argument('-bins', '--bins',
                        type=int, default=7,
                        help=help_str)
    help_str = 'Output file containing non-coding frameshift lengths'
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help=help_str)
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # calculate the number of non-coding bases for a single genome
    genome_len = get_genome_length(opts['genome_file'])
    black_list_len = get_black_list_length(opts['black_list'])
    non_coding_len = genome_len - black_list_len
    non_coding_len = int(opts['coverage'] * non_coding_len)  # adjust for coverage

    # read in indel mutations file and keep only samples that likely have
    # genuine non-coding indels (to prevent skewing of background rate)
    fs_df = pd.read_csv(opts['indels'], sep='\t')
    sample_indel_cts = fs_df['Tumor_Sample'].value_counts()
    non_coding_samples = sample_indel_cts[sample_indel_cts>=opts['threshold']]
    fs_df = fs_df[fs_df['Tumor_Sample'].isin(non_coding_samples.index)]

    # get count information for frameshift lengths
    fs_len, fs_count = get_frameshift_info(fs_df, bins=opts['bins'])

    # count number of distinct samples
    num_samples = len(fs_df['Tumor_Sample'].unique())

    # calculate total "bases at risk"
    bases_at_risk = non_coding_len * num_samples

    # format information
    tmp_info = [genome_len, black_list_len, non_coding_len,
                num_samples, bases_at_risk]
    out_info = fs_count + tmp_info

    # format header
    tmp_header = ['Genome Length', 'Black List Length', 'Non-coding Length',
                  'Number of Samples', 'Bases at Risk']
    tmp_count_header = map(str, fs_len)
    out_header = tmp_count_header + tmp_header

    # write to file
    out_df = pd.DataFrame([out_info],
                          columns=out_header,
                          index=['non-coding frameshift'])
    out_df.to_csv(opts['output'], sep='\t')



if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)

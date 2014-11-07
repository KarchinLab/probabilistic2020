import pandas as pd
import argparse


def correct_chrom_names(chroms):
    chrom_list = []
    for chrom in chroms:
        # fix chrom numbering
        chrom = chrom.replace('23', 'X')
        chrom = chrom.replace('24', 'Y')
        chrom = chrom.replace('25', 'Mt')
        if not chrom.startswith('chr'):
            chrom = 'chr' + chrom
        chrom_list.append(chrom)
    return chrom_list


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mutations',
                        required=True, type=str,
                        help='Unmapped mutations')
    parser.add_argument('-b', '--bed',
                        required=True, type=str,
                        help='Output bed file')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    df = pd.read_csv(opts['mutations'], sep='\t')
    df['Chromosome'] = correct_chrom_names(df['Chromosome'])
    df = df.sort(columns=['Chromosome', 'Start_Position', 'End_Position'])
    df[['Chromosome', 'Start_Position', 'End_Position']].to_csv(opts['bed'],
                                                                header=False,
                                                                sep='\t',
                                                                index=False)


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)

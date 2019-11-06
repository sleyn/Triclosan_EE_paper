import gff
import pandas as pd
import re
import argparse
import glob
from scipy.stats import chisquare
from statistics import mean

class timeline:
    def __init__(self, vcf_list, gff_file, ann_vcf, repeat_file):
        print('Initialize summary table')
        self.vcf_list = vcf_list
        self.gff_file = gff_file
        self.ann_vcf = ann_vcf
        self.repeat_file = repeat_file
        self.rregs = dict()  # repeat regions
        self.table = pd.DataFrame(columns=['Chrom',
                                           'Pos',
                                           'Ref',
                                           'Alt'])      # main output table
        self.ann_table = pd.DataFrame(columns=['Chrom',
                                               'Pos',
                                               'Ref',
                                               'Alt',
                                               'Effect',
                                               'Gene',
                                               'Locus_tag',
                                               'Mutation'])  # table with parsed snpEff output
        self.gff = gff.gff(self.gff_file)                    # create gff class object
        self.samples = self.collect_snames()                 # collect sample names

    # collect sample names
    def collect_snames(self):
          return [re.search('_([^_]+)\.vcf', path).group(1) for path in self.vcf_list]

    # process vcf files
    def process_vcf(self):
        print('Process VCF files for time points')
        for vcf_path in self.vcf_list:
            print('Process: ' + vcf_path)
            sname = re.search('_([^_]+)\.vcf', vcf_path).group(1) # exctract sample name
            with open(vcf_path, 'r') as vcf_file:
                temp_table = list() # temporary table for current sample
                for line in vcf_file.readlines():
                    if line[0] == '#':      # skip description lines
                        continue

                    var = line.split('\t')
                    # Chrom Pos Ref Alt Freq
                    temp_table.append([var[0], var[1], var[3], var[4], re.search('AF=([^;]+);', var[7]).group(1)])

                sample_table = pd.DataFrame(temp_table, columns=['Chrom',
                                                                 'Pos',
                                                                 'Ref',
                                                                 'Alt',
                                                                 sname])
                sample_table[sname] = sample_table[sname].astype(float)

                self.table = pd.merge(self.table, sample_table, how='outer', on=['Chrom', 'Pos', 'Ref', 'Alt'])   # merge temporary table to the main table

    # Process SnpEff annotation VCF file and populate table
    def snpEff(self):
        print('Add snpEff annotations')

        temp_table = list()

        with open(self.ann_vcf, 'r') as ann_file:
            for line in ann_file.readlines():
                if line[0] == '#':      # skip description lines
                    continue

                var = line.split('\t')
                eff = var[7]

                if 'ANN' not in eff:        # sometimes snpEff produces no annotations
                    print("!!!!!!Missing annotation!!!!!!")
                    print(line)
                    continue

                eff = re.sub('.+;ANN', 'ANN', eff)
                eff = re.sub('fig\|\d+', '', eff)
                effs = eff.split('|')

                # Chrom Pos Ref Alt Effect Gene Locus_tag Mutation
                temp_table.append([var[0],
                                   var[1],
                                   var[3],
                                   var[4],
                                   effs[1],
                                   effs[3],
                                   effs[4],
                                   effs[10][2:]])

            ann_table = pd.DataFrame(temp_table, columns=['Chrom',
                                               'Pos',
                                               'Ref',
                                               'Alt',
                                               'Effect',
                                               'Gene',
                                               'Locus_tag',
                                               'Mutation'])

            self.ann_table = self.ann_table.append(ann_table)

        self.table = pd.merge(self.table, self.ann_table, how='left', on=['Chrom', 'Pos', 'Ref', 'Alt']) # merge annotations from snpEff to the main table

    # add product annotation from gff file
    def add_ann(self):
        print('Add product annotations')
        for i in range(len(self.table)):
            self.table.at[i, 'Annotation'] = self.gff.gff_pos[self.table.iloc[i]['Chrom']][int(self.table.iloc[i]['Pos'])][3]

    # Auxillary chi^2 test function
    @staticmethod
    def _chisq(x):
        x = [v*100 for v in x]      # transform to percents
        if mean(x) == 0:
            x = [1 for value in x]

        return chisquare(x)

    # make ch.sq test for uniform distribution
    def test_uniform(self):
        # find how many reactors do we have
        print('Make Chi^2 uniformity test')
        reactors = sorted(set([re.search('(\d+)\w', sample).group(1) for sample in self.samples]))

        for rtr in reactors:
            tps = [tp for tp in self.samples if re.search(rtr, tp)]        # collect time points for reactor
            self.table['Uniform_' + rtr] = self.table[tps].apply(lambda x: self._chisq(x.tolist())[1], axis=1)   # Chi square test for uniformity

    # read repeat regions file:
    def _read_rep(self):
        with open(self.repeat_file, 'r') as r_file:
            for rep in r_file.readlines():
                rep = rep.rstrip().split('\t')
                if rep[0] in self.rregs:
                    self.rregs[rep[0]].append([int(rep[1]), int(rep[2])])  # chrom => [start, end]
                else:
                    self.rregs[rep[0]] = [[int(rep[1]), int(rep[2])]]

    # test id current variant is in repeat
    def _test_repvar(self, cp):
        for reg in self.rregs[cp['Chrom']]:
            if reg[0] <= int(cp['Pos']) <= reg[1]:
                return 'IN_REPEAT'
        return 'NOT'

    # test if variant is in a repeat region
    def test_repeat(self):
        print('Test if variants are in repeat regions')
        self._read_rep()
        self.table['IN_REPEAT'] = self.table[['Chrom', 'Pos']].apply(lambda x: self._test_repvar(x), axis=1)


# executable part

parser = argparse.ArgumentParser(description='Script combines VCF outputs for time series in summary table')
parser.add_argument('-g', '--gff', action='store', type=str, help='GFF file for reference')
parser.add_argument('-a', '--ann', action='store', type=str, help='merged snpEff VCF file')
parser.add_argument('-r', '--rep', action='store', type=str, help='file with repeat regions')
parser.add_argument('-v', '--vcf_folder', action='store', type=str, help='folder with vcf files from time series samples')
parser.add_argument('-o', '--out', action='store', type=str, help='output file')
args = parser.parse_args()

vcf_path = args.vcf_folder
print(args.vcf_folder)
if vcf_path[-1] != '/':     # add final '/' if there is none
    vcf_path += '/'
vcf_files = glob.glob(vcf_path + '*[0-9][A-Z].vcf')

table = timeline(vcf_files, args.gff, args.ann, args.rep)
table.process_vcf()
table.snpEff()
print('Process GFF file: ' + args.gff)
table.gff.readgff()
table.add_ann()
table.table = table.table.fillna(0)
table.test_uniform()
table.test_repeat()
table.table.to_csv(args.out, sep='\t', index=False, quoting=3)

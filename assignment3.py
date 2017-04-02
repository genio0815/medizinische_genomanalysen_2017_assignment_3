#! /usr/bin/env python2

import vcf
import vcf.utils
import hgvs
import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.assemblymapper
import os
from bioutils.assemblies import make_name_ac_map
import logging

__author__ = 'alexander.bindeus'


class Assignment3:
    
    def __init__(self):
        print("PyVCF version: %s" % vcf.VERSION)
        print("HGVS version: %s" % hgvs.__version__)

        # hard coded input...do not care
        self.sonFile = os.path.join(os.getcwd(), "AmpliseqExome.20141120.NA24385.vcf")
        self.motherFile = os.path.join(os.getcwd(), "AmpliseqExome.20141120.NA24143.vcf")
        self.fatherFile = os.path.join(os.getcwd(), "AmpliseqExome.20141120.NA24149.vcf")

    def get_total_number_of_variants(self, printName, infile):
        vcf_reader = vcf.Reader(open(infile, 'rb'))
        count = 0
        for record in vcf_reader:
            count += 1

        print "variants in", printName, "file:\t", count

    def get_variants_shared_by_two(self, name1, name2, file1, file2):
        readerA = vcf.Reader(open(file1, 'rb'))
        readerB = vcf.Reader(open(file2, 'rb'))

        count = 0

        for a, b in vcf.utils.walk_together(readerA, readerB):
            if a and b \
                    and a.ALT == b.ALT \
                    and a.CHROM == b.CHROM \
                    and a.INFO.get("GT") == b.INFO.get("GT") \
                    and a.POS == b.POS \
                    and a.REF == b.REF:
                count += 1

        print "variants shared by", name1, "and", name2, ":\t", count

    def get_variants_shared_by_trio(self):
        readerA = vcf.Reader(open(self.sonFile, 'rb'))
        readerB = vcf.Reader(open(self.motherFile, 'rb'))
        readerC = vcf.Reader(open(self.fatherFile, 'rb'))

        count = 0

        for line in vcf.utils.walk_together(readerA, readerB, readerC):
            a = line[0]
            b = line[1]
            c = line[2]

            if a and b and c \
                    and a.ALT == b.ALT == c.ALT \
                    and a.CHROM == b.CHROM == c.CHROM \
                    and a.INFO.get("GT") == b.INFO.get("GT") == c.INFO.get("GT") \
                    and a.POS == b.POS == c.POS \
                    and a.REF == b.REF == c.REF:
                count += 1

        print "variants shared by father, mother and son", count
        
    def merge_mother_father_son_into_one_vcf(self):
        readerA = vcf.Reader(open(self.sonFile, 'rb'))
        readerB = vcf.Reader(open(self.motherFile, 'rb'))
        readerC = vcf.Reader(open(self.fatherFile, 'rb'))

        f = open(os.path.join(os.getcwd(), 'mergedVariants.vcf'), 'w')
        out = vcf.Writer(f, readerA, '\n')

        for line in vcf.utils.walk_together(readerA, readerB, readerC):
            for entry in line:
                if entry != None:  # no duplicates
                    out.write_record(entry)
                    break

        out.close()

    def convert_first_variants_of_son_into_HGVS(self):

        print("converting first 100 varaints of son to HGFS format")

        ## Connect to UTA
        hdp = hgvs.dataproviders.uta.connect()
        logging.basicConfig()
        assembly_mapper = hgvs.assemblymapper.AssemblyMapper(hdp, normalize=False)  # EasyVariantMapper before
        ## Used for parsing
        hgvsparser = hgvs.parser.Parser()  # Parser

        reader = vcf.Reader(open(self.sonFile, 'rb'))
        outfile = open("first_100_variants_son.out", "w")

        def mapping(genome_hgvs):
            g = hgvsparser.parse_hgvs_variant(genome_hgvs)
            for tr in assembly_mapper.relevant_transcripts(g):
                try:
                    c = assembly_mapper.g_to_c(g, tr)  # coding
                    outfile.writelines("%s\t%s\n" % (g, c))
                except hgvs.exceptions.HGVSUsageError:
                    n = assembly_mapper.g_to_n(g, tr)  # non coding
                    outfile.writelines("%s\t%s\n" % (g, n))
                except hgvs.exceptions.HGVSInvalidIntervalError:
                    outfile.writelines("mapping error at %s\t%s\n" % (g, tr))

        limit = 100
        count = 0

        for record in reader:
            if count < limit:
                refseq_nc_number = make_name_ac_map("GRCh37.p13")[record.CHROM[3:]]
                try:
                    genome_hgvs = "%s:g.%s%s>%s" % (refseq_nc_number, str(record.POS), str(record.REF),
                                                    str(record.ALT[0]))
                    mapping(genome_hgvs)
                except Exception as e:
                    print("caught exception", e)
            else:
                break

            count += 1

        print("Wrote first 100 variants of son file into file '"' first_100_variants_son.out'"' ")
        outfile.close()
    
    def print_summary(self):
        self.get_total_number_of_variants('mother', self.motherFile)
        self.get_total_number_of_variants('father', self.fatherFile)
        self.get_total_number_of_variants('son', self.sonFile)

        self.get_variants_shared_by_two('mother', 'son', self.motherFile, self.sonFile)
        self.get_variants_shared_by_two('father', 'son', self.fatherFile, self.sonFile)
        self.get_variants_shared_by_trio()
        self.convert_first_variants_of_son_into_HGVS()
        self.merge_mother_father_son_into_one_vcf()


if __name__ == '__main__':
    print("Assignment 3")
    print(__author__)
    assignment1 = Assignment3()
    assignment1.print_summary()
    
    


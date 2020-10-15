#!/usr/bin/env python

# Standard library imports
import os
import sys
import re
import logging

# Third party imports
import argparse
import subprocess
import datetime


# Local application imports
from misc import check_file_exists, extract_sample, check_create_dir, execute_subprocess, \
    extract_read_list, file_to_list, obtain_group_cov_stats, clean_unwanted_files, \
    check_reanalysis, vcf_stats, remove_low_quality, obtain_overal_stats
from preprocessing import fastqc_quality, fastp_trimming, format_html_image
from pe_mapper import bwa_mapping, sam_to_index_bam
from bam_variant import picard_dictionary, samtools_faidx, picard_markdup, ivar_trim, ivar_variants, ivar_consensus, \
    replace_consensus_header, create_bamstat, create_coverage
from vcf_process import filter_tsv_variants, vcf_consensus_filter, highly_hetz_to_bed, poorly_covered_to_bed, non_genotyped_to_bed
from annotation import annotate_snpeff, annotate_pangolin
from species_determination import mash_screen, extract_species_from_screen
from compare_snp import ddtb_add, ddtb_compare, recalibrate_ddbb_vcf

"""
=============================================================
HEADER
=============================================================

INSTITUTION:IiSGM
AUTHOR: Pedro J. Sola (pedroscampoy@gmail.com)
d^v^b
VERSION=0.1
CREATED: 22 Sep 2020
REVISION: 


TODO:
    Adapt check_reanalysis
    Check file with multiple arguments
    Check program is installed (dependencies)
================================================================
END_OF_HEADER
================================================================
"""

#COLORS AND AND FORMATTING

END_FORMATTING = '\033[0m'
WHITE_BG = '\033[0;30;47m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
RED = '\033[31m'
GREEN = '\033[32m'
MAGENTA = '\033[35m'
BLUE =  '\033[34m'
CYAN = '\033[36m'
YELLOW = '\033[93m'
DIM = '\033[2m'

logger = logging.getLogger()

def main():
    """
    Create main function to capture code errors: https://stackoverflow.com/questions/6234405/logging-uncaught-exceptions-in-python
    """

    #ARGUMENTS

    def get_arguments():

        parser = argparse.ArgumentParser(prog = 'covidma.py', description= 'Pipeline to call variants (SNVs) with any non model organism. Specialised in SARS-CoV-2')
        
        input_group = parser.add_argument_group('Input', 'Input parameters')

        input_group.add_argument('-i', '--input', dest="input_dir", metavar="input_directory", type=str, required=True, help='REQUIRED.Input directory containing all fast[aq] files')
        input_group.add_argument('-r', '--reference', metavar="reference", type=str, required=True, help='REQUIRED. File to map against')
        input_group.add_argument('-a', '--annotation', metavar="annotation", type=str, required=True, help='REQUIRED. gff3 file to annotate variants')
        input_group.add_argument('-s', '--sample', metavar="sample", type=str, required=False, help='Sample to identify further files')
        input_group.add_argument('-S', '--sample_list', type=str, required=False, help='Sample names to analyse only in the file supplied')
        input_group.add_argument('-p', '--primers', type=str, default='/home/laura/COVID/primers/nCoV-2019.bed', required=False, help='Bed file including primers to trim')
        input_group.add_argument('-B', '--annot_bed', type=str, required=False, action='append', help='bed file to annotate')
        input_group.add_argument('-V', '--annot_vcf', type=str, required=False, action='append', help='vcf file to annotate')

        quality_group = parser.add_argument_group('Quality parameters', 'parameters for diferent triming conditions')

        quality_group.add_argument('-c', '--coverage20', type=int, default=90, required=False, help='Minimum percentage of coverage at 20x to clasify as uncovered (Default 90)')
        quality_group.add_argument('-n', '--min_snp', type=int, required=False, default=1, help='SNP number to pass quality threshold')

        output_group = parser.add_argument_group('Output', 'Required parameter to output results')

        output_group.add_argument('-o', '--output', type=str, required=True, help='REQUIRED. Output directory to extract all results')
        output_group.add_argument('-C', '--noclean', required=False, action='store_false', help='Clean unwanted files for standard execution')

        params_group = parser.add_argument_group('Parameters', 'parameters for diferent stringent conditions')

        params_group.add_argument('-T', '--threads', type=str, dest = "threads", required=False, default=16, help='Threads to use')
        params_group.add_argument('-M', '--memory', type=str, dest = "memory", required=False, default=32, help='Max memory to use')
        
        """
        vcf_group = parser.add_argument_group('VCF filters', 'parameters for variant filtering')

        vcf_group.add_argument('-b', '--bed_remove', type=str, required=False, default=False, help='BED file with position ranges to filter from final vcf')
        vcf_group.add_argument('-m', '--maxnocallfr', type=str, required=False, default=0.1, help='maximun proportion of samples with non genotyped alleles')

        annot_group = parser.add_argument_group('Annotation', 'parameters for variant annotation')

        annot_group.add_argument('--mash_database', type=str, required=False, default="/home/pjsola/REFERENCES/mash/RefSeq88n.msh", help='MASH ncbi annotation containing all species database')
        annot_group.add_argument('--snpeff_database', type=str, required=False, default=False, help='snpEFF annotation database')
        """
        
    

        arguments = parser.parse_args()

        return arguments

    args = get_arguments()


    ######################################################################
    #####################START PIPELINE###################################
    ######################################################################
    output = os.path.abspath(args.output)
    group_name = output.split("/")[-1]
    reference = os.path.abspath(args.reference)
    annotation = os.path.abspath(args.annotation)

    #LOGGING
    #Create log file with date and time
    right_now = str(datetime.datetime.now())
    right_now_full = "_".join(right_now.split(" "))
    log_filename = group_name + "_" + right_now_full + ".log"
    log_folder = os.path.join(output, 'Logs')
    check_create_dir(log_folder)
    log_full_path = os.path.join(log_folder, log_filename)

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s:%(message)s')

    file_handler = logging.FileHandler(log_full_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    #stream_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)


    logger.info("\n\n" + BLUE + BOLD + "STARTING PIPELINE IN GROUP: " + group_name + END_FORMATTING)

    today = str(datetime.date.today())

    logger.info("ARGUMENTS:")
    logger.info(str(args))

    check_reanalysis(args.output)

    #Obtain all R1 and R2 from folder
    r1, r2 = extract_read_list(args.input_dir)

    #Check if there are samples to filter out
    sample_list_F = []
    if args.sample_list == None:
        logger.info("\n" + "No samples to filter")
        for r1_file, r2_file in zip(r1, r2):
            sample = extract_sample(r1_file, r2_file)
            sample_list_F.append(sample)
    else:
        logger.info("samples will be filtered")
        sample_list_F = file_to_list(args.sample_list)
    logger.info("\n%d samples will be analysed: %s" % (len(sample_list_F), ",".join(sample_list_F)))


    #PREPARE REFERENCE FOR MAPPING + FAI + DICT #########
    #####################################################
    
    #picard_dictionary(args)
    samtools_faidx(args)

    #DECLARE FOLDERS CREATED IN PIPELINE ################
    #AND KEY FILES ######################################
    #####################################################
    #Annotation related parameters
    #script_dir = os.path.dirname(os.path.realpath(__file__))

    #Output related
    out_qc_dir = os.path.join(output, "Quality")
    out_qc_pre_dir = os.path.join(out_qc_dir, "raw") #subfolder
    out_qc_post_dir = os.path.join(out_qc_dir, "processed") #subfolder
    out_trim_dir = os.path.join(output, "Trimmed")
    out_map_dir = os.path.join(output, "Bam")
    out_variant_dir = os.path.join(output, "Variants")
    out_variant_ivar_dir = os.path.join(out_variant_dir, "ivar_raw") #subfolder
    out_filtered_ivar_dir = os.path.join(out_variant_dir, "ivar_filtered") #subfolder
    out_consensus_dir = os.path.join(output, "Consensus")

    out_stats_dir = os.path.join(output, "Stats")
    out_stats_bamstats_dir = os.path.join(out_stats_dir, "Bamstats") #subfolder
    out_stats_coverage_dir = os.path.join(out_stats_dir, "Coverage") #subfolder
    out_compare_dir = os.path.join(output, "Compare")

    out_annot_dir = os.path.join(output, "Annotation")
    out_annot_snpeff_dir = os.path.join(out_annot_dir, "snpeff")
    out_annot_pangolin_dir = os.path.join(out_annot_dir, "pangolin")

    #highly_hetz_bed = os.path.join(out_variant_dir, "highly_hetz.bed")
    #non_genotyped_bed = os.path.join(out_variant_dir, "non_genotyped.bed")
    #poorly_covered_bed = os.path.join(out_cov_dir, "poorly_covered.bed")

    for r1_file, r2_file in zip(r1, r2):
        #EXtract sample name
        sample = extract_sample(r1_file, r2_file)
        args.sample = sample
        if sample in sample_list_F:

            sample_number = str(sample_list_F.index(sample) + 1)
            sample_total = str(len(sample_list_F))

            out_markdup_trimmed_name = sample + ".rg.markdup.trimmed.sorted.bam"
            output_markdup_trimmed_file = os.path.join(out_map_dir, out_markdup_trimmed_name)

            logger.info("\n" + WHITE_BG + "STARTING SAMPLE: " + sample + " (" + sample_number + "/" + sample_total + ")" + END_FORMATTING)

            if not os.path.isfile(output_markdup_trimmed_file):
            
                args.r1_file = r1_file
                args.r2_file = r2_file

                ##############START PIPELINE#####################
                #################################################

                #INPUT ARGUMENTS
                ################
                check_file_exists(r1_file)
                check_file_exists(r2_file)

                args.output = os.path.abspath(args.output)
                check_create_dir(args.output)

                ######################QUALITY CHECK in RAW with fastqc
                ######################################################
                check_create_dir(out_qc_dir)

                out_qc_raw_name_r1 = (".").join(r1_file.split('/')[-1].split('.')[0:-2]) + '_fastqc.html'
                out_qc_raw_name_r2 = (".").join(r2_file.split('/')[-1].split('.')[0:-2]) + '_fastqc.html'
                output_qc_raw_file_r1 = os.path.join(out_qc_pre_dir, out_qc_raw_name_r1)
                output_qc_raw_file_r2 = os.path.join(out_qc_pre_dir, out_qc_raw_name_r2)
                                
                if os.path.isfile(output_qc_raw_file_r1) and os.path.isfile(output_qc_raw_file_r2):
                    logger.info(YELLOW + DIM + output_qc_raw_file_r1 + " EXIST\nOmmiting QC for sample " + sample + END_FORMATTING)
                else:
                    logger.info(GREEN + "Checking quality in sample " + sample + END_FORMATTING)
                    logger.info("R1: " + r1_file + "\nR2: " + r2_file)
                    fastqc_quality(r1_file, r2_file, out_qc_pre_dir, args.threads)

                """
                TODO: Human filter
                """
                
                ####QUALITY TRIMMING AND ADAPTER REMOVAL WITH fastp
                ###################################################
                out_trim_name_r1 = sample + ".trimmed_R1.fastq.gz"
                out_trim_name_r2 = sample + ".trimmed_R2.fastq.gz"
                output_trimming_file_r1 = os.path.join(out_trim_dir, out_trim_name_r1)
                output_trimming_file_r2 = os.path.join(out_trim_dir, out_trim_name_r2)
                
                if os.path.isfile(output_trimming_file_r1) and os.path.isfile(output_trimming_file_r2):
                    logger.info(YELLOW + DIM + output_trimming_file_r1 + " EXIST\nOmmiting Trimming for sample " + sample + END_FORMATTING)
                else:
                    logger.info(GREEN + "Trimming sample " + sample + END_FORMATTING)
                    fastp_trimming(r1_file, r2_file, sample, out_trim_dir, threads=args.threads, min_qual=20, window_size=10, min_len=35)


                ##################QUALITY CHECK in TRIMMED with fastqc
                ######################################################
                check_create_dir(out_qc_dir)
                
                out_qc_pos_r1 = sample + ".trimmed_R1_fastqc.html"
                out_qc_pos_r2 = sample + ".trimmed_R2_fastqc.html"
                output_qc_precessed_file_r1 = os.path.join(out_qc_post_dir, out_qc_pos_r1)
                output_qc_precessed_file_r2 = os.path.join(out_qc_post_dir, out_qc_pos_r2)
                                
                if os.path.isfile(output_qc_precessed_file_r1) and os.path.isfile(output_qc_precessed_file_r2):
                    logger.info(YELLOW + DIM + output_qc_raw_file_r1 + " EXIST\nOmmiting QC for sample " + sample + END_FORMATTING)
                else:
                    logger.info(GREEN + "Checking quality in processed sample " + sample + END_FORMATTING)
                    logger.info("R1: " + output_trimming_file_r1 + "\nR2: " + output_trimming_file_r2)
                    fastqc_quality(output_trimming_file_r1, output_trimming_file_r2, out_qc_post_dir, args.threads)

                #MAPPING WITH BWA - SAM TO SORTED BAM - ADD HEADER SG
                #####################################################
                out_map_name = sample + ".rg.sorted.bam"
                output_map_file = os.path.join(out_map_dir, out_map_name)

                if os.path.isfile(output_map_file):
                    logger.info(YELLOW + DIM + output_map_file + " EXIST\nOmmiting Mapping for sample " + sample + END_FORMATTING)
                else:
                    logger.info(GREEN + "Mapping sample " + sample + END_FORMATTING)
                    logger.info("R1: " + output_trimming_file_r1 + "\nR2: " + output_trimming_file_r2 + "\nReference: " + reference)
                    bwa_mapping(output_trimming_file_r1, output_trimming_file_r2, reference, sample, out_map_dir, threads=args.threads)
                    sam_to_index_bam(sample, out_map_dir, output_trimming_file_r1, threads=args.threads)


                #MARK DUPLICATES WITH PICARDTOOLS ###################
                #####################################################
                out_markdup_name = sample + ".rg.markdup.sorted.bam"
                output_markdup_file = os.path.join(out_map_dir, out_markdup_name)

                if os.path.isfile(output_markdup_file):
                    logger.info(YELLOW + DIM + output_markdup_file + " EXIST\nOmmiting Duplucate Mark for sample " + sample + END_FORMATTING)
                else:
                    logger.info(GREEN + "Marking Dupes in sample " + sample + END_FORMATTING)
                    logger.info("Input Bam: " + output_map_file)
                    picard_markdup(output_map_file)
                
                #TRIM PRIMERS WITH ivar trim ########################
                #####################################################

                if os.path.isfile(output_markdup_trimmed_file):
                    logger.info(YELLOW + DIM + output_markdup_trimmed_file + " EXIST\nOmmiting Duplucate Mark for sample " + sample + END_FORMATTING)
                else:
                    logger.info(GREEN + "Trimming primers in sample " + sample + END_FORMATTING)
                    logger.info("Input Bam: " + output_markdup_file)
                    ivar_trim(output_markdup_file, args.primers, sample, min_length=30, min_quality=20, sliding_window_width=4)
            else:
                logger.info(YELLOW + DIM + output_markdup_trimmed_file + " EXIST\nOmmiting BAM mapping and BAM manipulation in sample " + sample + END_FORMATTING)
            
            ########################END OF MAPPING AND BAM MANIPULATION#####################################################################
            ################################################################################################################################
            
            #VARIANT CALLING WTIH ivar variants##################
            #####################################################
            check_create_dir(out_variant_dir)
            out_ivar_variant_name = sample + ".tsv"
            out_ivar_variant_file = os.path.join(out_variant_ivar_dir, out_ivar_variant_name)

            if os.path.isfile(out_ivar_variant_file):
                logger.info(YELLOW + DIM + out_ivar_variant_file + " EXIST\nOmmiting Variant call for  sample " + sample + END_FORMATTING)
            else:
                logger.info(GREEN + "Calling variants with ivar in sample " + sample + END_FORMATTING)
                ivar_variants(reference, output_markdup_trimmed_file, out_variant_dir, sample, annotation, min_quality=20, min_frequency_threshold=0.05, min_depth=5)


            #VARIANT FILTERING ##################################
            #####################################################
            check_create_dir(out_filtered_ivar_dir)
            out_ivar_filtered_file = os.path.join(out_filtered_ivar_dir, out_ivar_variant_name)

            if os.path.isfile(out_ivar_filtered_file):
                logger.info(YELLOW + DIM + out_ivar_filtered_file + " EXIST\nOmmiting Variant filtering for  sample " + sample + END_FORMATTING)
            else:
                logger.info(GREEN + "Filtering variants in sample " + sample + END_FORMATTING)
                filter_tsv_variants(out_ivar_variant_file, out_filtered_ivar_dir, min_frequency=0.8, min_total_depth=10, min_alt_dp=4, is_pass=True, only_snp=True)
            
            #CREATE CONSENSUS with ivar consensus##################
            #######################################################
            check_create_dir(out_consensus_dir)
            out_ivar_consensus_name = sample + ".fa"
            out_ivar_consensus_file = os.path.join(out_consensus_dir, out_ivar_consensus_name)

            if os.path.isfile(out_ivar_consensus_file):
                logger.info(YELLOW + DIM + out_ivar_consensus_file + " EXIST\nOmmiting Consensus for  sample " + sample + END_FORMATTING)
            else:
                logger.info(GREEN + "Creating consensus with ivar in sample " + sample + END_FORMATTING)
                ivar_consensus(output_markdup_trimmed_file, out_consensus_dir, sample, min_quality=20, min_frequency_threshold=0.8, min_depth=20, uncovered_character='N')
                logger.info(GREEN + "Replacing consensus header in " + sample + END_FORMATTING)
                replace_consensus_header(out_ivar_consensus_file)


            ########################CREATE STATS AND QUALITY FILTERS########################################################################
            ################################################################################################################################
            #CREATE Bamstats#######################################
            #######################################################
            check_create_dir(out_stats_dir)
            check_create_dir(out_stats_bamstats_dir)
            out_bamstats_name = sample + ".bamstats"
            out_bamstats_file = os.path.join(out_stats_bamstats_dir, out_bamstats_name)

            if os.path.isfile(out_bamstats_file):
                logger.info(YELLOW + DIM + out_bamstats_file + " EXIST\nOmmiting Bamstats for  sample " + sample + END_FORMATTING)
            else:
                logger.info(GREEN + "Creating bamstats in sample " + sample + END_FORMATTING)
                create_bamstat(output_markdup_trimmed_file, out_stats_bamstats_dir, sample, threads=args.threads)

            #CREATE Bamstats#######################################
            #######################################################
            check_create_dir(out_stats_coverage_dir)
            out_coverage_name = sample + ".cov"
            out_coverage_file = os.path.join(out_stats_coverage_dir, out_coverage_name)

            if os.path.isfile(out_coverage_file):
                logger.info(YELLOW + DIM + out_coverage_file + " EXIST\nOmmiting Bamstats for  sample " + sample + END_FORMATTING)
            else:
                logger.info(GREEN + "Creating coverage in sample " + sample + END_FORMATTING)
                create_coverage(output_markdup_trimmed_file, out_stats_coverage_dir, sample)

            
    
    ###################fastqc OUTPUT FORMAT FOR COMPARISON
    ######################################################
    logger.info(GREEN + "Creating summary report for quality result " + END_FORMATTING)
    format_html_image(out_qc_dir)

    ###############################coverage OUTPUT SUMMARY
    ######################################################
    logger.info(GREEN + "Creating summary report for coverage result " + END_FORMATTING)
    obtain_group_cov_stats(out_stats_coverage_dir, group_name)

    #####################READS and VARIANTS OUTPUT SUMMARY
    ######################################################
    logger.info(GREEN + "Creating overal summary report " + END_FORMATTING)
    obtain_overal_stats(output, group_name)

    ######################################REMOVE UNCOVERED
    ##############################################################################################################################
    logger.info(GREEN + "Removing low quality samples" + END_FORMATTING)
    remove_low_quality(output, min_percentage_20x=args.coverage20, min_hq_snp=args.min_snp, type_remove='Uncovered')

    #ANNOTATION WITH SNPEFF aND PANGOLIN ################
    #####################################################
    logger.info("\n\n" + BLUE + BOLD + "STARTING ANNOTATION IN GROUP: " + group_name + END_FORMATTING + "\n")
    check_create_dir(out_annot_dir)
    check_create_dir(out_annot_snpeff_dir)
    check_create_dir(out_annot_pangolin_dir)
    ####SNPEFF
    for root, _, files in os.walk(out_filtered_ivar_dir):
        if root == out_filtered_ivar_dir: 
            for name in files:
                if name.endswith('.tsv'):
                    sample = name.split('.')[0]
                    filename = os.path.join(root, name)
                    out_annot_file = os.path.join(out_annot_snpeff_dir, sample + ".annot")
                    if os.path.isfile(out_annot_file):
                        logger.info(YELLOW + DIM + out_annot_file + " EXIST\nOmmiting snpEff Annotation for sample " + sample + END_FORMATTING)
                    else:
                        logger.info(GREEN + "Annotating sample with snpEff: " + sample + END_FORMATTING)
                        output_vcf = os.path.join(out_annot_snpeff_dir, sample + '.vcf')
                        annotate_snpeff(filename, output_vcf, out_annot_file)
    
    ####PANGOLIN
    for root, _, files in os.walk(out_consensus_dir):
        if root == out_consensus_dir: 
            for name in files:
                if name.endswith('.fa'):
                    sample = name.split('.')[0]
                    filename = os.path.join(root, name)
                    out_pangolin_filename = sample + ".lineage.csv"
                    out_pangolin_file = os.path.join(out_annot_pangolin_dir, out_pangolin_filename)
                    if os.path.isfile(out_pangolin_file):
                        logger.info(YELLOW + DIM + out_pangolin_file + " EXIST\nOmmiting Lineage for  sample " + sample + END_FORMATTING)
                    else:
                        logger.info(GREEN + "Obtaining Lineage in sample " + sample + END_FORMATTING)
                        annotate_pangolin(filename, out_annot_pangolin_dir, out_pangolin_filename, threads=args.threads, max_ambig=0.6)









    ################SNP COMPARISON using tsv variant files
    ######################################################
    logger.info("\n\n" + BLUE + BOLD + "STARTING COMPARISON IN GROUP: " + group_name + END_FORMATTING + "\n")

    check_create_dir(out_compare_dir)
    folder_compare = today + "_" + group_name
    path_compare = os.path.join(out_compare_dir, folder_compare)
    check_create_dir(path_compare)
    full_path_compare = os.path.join(path_compare, group_name)
    compare_snp_matrix = full_path_compare + ".tsv"

    ddtb_add(out_filtered_ivar_dir, full_path_compare)
    compare_snp_matrix_recal = full_path_compare + ".revised.tsv"
    recalibrated_snp_matrix = recalibrate_ddbb_vcf(compare_snp_matrix, out_map_dir)
    recalibrated_snp_matrix.to_csv(compare_snp_matrix_recal, sep="\t", index=False)
    ddtb_compare(compare_snp_matrix_recal, distance=0)

    logger.info("\n\n" + MAGENTA + BOLD + "COMPARING FINISHED IN GROUP: " + group_name + END_FORMATTING + "\n")

    """
                
                
                
                
                #CALCULATE COVERAGE FOR EACH POSITION##################
                #######################################################
                out_cov_name = sample + ".cov"
                output_cov_file = os.path.join(out_cov_dir, out_cov_name)

                if os.path.isfile(output_cov_file):
                    logger.info(YELLOW + DIM + output_cov_file + " EXIST\nOmmiting coverage calculation for sample " + sample + END_FORMATTING)
                else:
                    logger.info(GREEN + "Calculating coverage in sample " + sample + END_FORMATTING)
                    get_coverage(args, output_markdup_file, output_fmt="-d")

                #SPECIES DETERMINATION USING mash #################
                ###################################################
                out_mash_name = sample + ".screen.tab"
                output_mash_file = os.path.join(out_species_dir, out_mash_name)
                
                if os.path.isfile(output_mash_file):
                    logger.info(YELLOW + DIM + output_mash_file + " EXIST\nOmmiting Species calculation for sample " + sample + END_FORMATTING)
                else:
                    logger.info(GREEN + "Determining species content in sample " + sample + END_FORMATTING)
                    mash_screen(args, winner=True, r2=False, mash_database=args.mash_database)
                

                #HAPLOTYPE CALL 1/2 FOR HARD FILTERING AND RECALIBRATION
                #######################################################
                out_gvcfr_name = sample + ".g.vcf"
                output_gvcfr_file = os.path.join(out_gvcfr_dir, out_gvcfr_name)

                if os.path.isfile(output_gvcfr_file):
                    logger.info(YELLOW + DIM + output_gvcfr_file + " EXIST\nOmmiting Haplotype Call (Recall) for sample " + sample + END_FORMATTING)
                else:
                    logger.info(GREEN + "Haplotype Calling (Recall) in sample " + sample + END_FORMATTING)
                    haplotype_caller(args, recalibrate=True, ploidy=2, bamout=False, forceactive=False)
            
            else:
                logger.info(YELLOW + DIM + "\nOMMITING BAM HANDLING FOR SAMPLE " + sample + END_FORMATTING)



    #GROUP COVERAGE SUMMARY STATS##########################
    #######################################################
    group_name = output.split("/")[-1]
    logger.info("\n\n" + BLUE + BOLD + "CHECKING LOW COVERED SAMPLES IN GROUP: " + group_name + END_FORMATTING + "\n")

    out_cov_name = group_name + ".coverage.tab"
    output_cov_file = os.path.join(out_cov_dir, out_cov_name)

    if os.path.isfile(output_cov_file):
        logger.info(YELLOW + DIM + output_cov_file + " EXIST\nOmmiting group coverage calculation for group " + group_name + END_FORMATTING)
        saples_low_covered = []
    else:
        logger.info(GREEN + "Group coverage stats in group " + group_name + END_FORMATTING)
        saples_low_covered = obtain_group_cov_stats(out_cov_dir, low_cov_threshold=args.mincov, unmmaped_threshold=args.unmmaped)


    if os.path.isfile(poorly_covered_bed):
        logger.info(YELLOW + DIM + output_cov_file + " EXIST\nOmmiting poorly covered calculation for group " + group_name + END_FORMATTING)
    else:
        logger.info(GREEN + "Calculating low covered regions " + group_name + END_FORMATTING)
        poorly_covered_to_bed(out_cov_dir, "poorly_covered", reference="CHROM", min_coverage=2, nocall_fr=0.5)

    if len(saples_low_covered) > 0:
        logger.info("\n" + YELLOW + BOLD + "There are sample(s) with low coverage that will be removed from the analysis: " + "\n"\
            + ",".join(saples_low_covered) + END_FORMATTING + "\n")
        remove_low_covered_mixed(args.output, saples_low_covered, "Uncovered")
        #Remove sample from the list of filtered samples
        ################################################
        for samples_to_remove in saples_low_covered:
            sample_list_F.remove(samples_to_remove)
    else:
        logger.info("\n" + YELLOW + BOLD + "All samples have a decent depth of coverage according to threshold supplied" + "\n")




    #ONCE ALL GVCF VARIANTS ARE CALLED, THEY ARE GATHERED AND FILTERED 
    # TO RECALIBRATE ORIGINAL MARKDUPPED BAM
    ######################################################################
    ##############START GROUP CALLING FOR RECALIBRATION###################
    ######################################################################

    group_name = output.split("/")[-1]
    logger.info("\n\n" + BLUE + BOLD + "STARTING JOINT CALL FOR RECALIBATION IN GROUP: " + group_name + END_FORMATTING + "\n")

    #CALL VARIANTS 1/2 FOR HARD FILTERING AND RECALIBRATION
    #######################################################
    out_gvcfr_name = group_name + ".cohort.g.vcf"
    output_gvcfr_file = os.path.join(out_gvcfr_dir, out_gvcfr_name)

    if os.path.isfile(output_gvcfr_file):
        logger.info(YELLOW + DIM + output_gvcfr_file + " EXIST\nOmmiting GVCF Combination (Recall) for group " + group_name + END_FORMATTING)
    else:
        logger.info(GREEN + "GVCF Combination (Recall) in group " + group_name + END_FORMATTING)
        combine_gvcf(args, recalibrate=True, all_gvcf=False)

    #CALL VARIANTS 1/2 FOR HARD FILTERING AND RECALIBRATION
    #######################################################
    out_vcfr_name = group_name + ".cohort.raw.vcf"
    output_vcfr_file = os.path.join(out_vcfr_dir, out_vcfr_name)

    if os.path.isfile(output_vcfr_file):
        logger.info(YELLOW + DIM + output_vcfr_file + " EXIST\nOmmiting Variant Calling (Recall-Group) for group " + group_name + END_FORMATTING)
    else:
        logger.info(GREEN + "Variant Calling (Recall-Group) in group " + group_name + END_FORMATTING)
        call_variants(args, recalibrate=True, group=True)

    #SELECT VARIANTS 1/2 FOR HARD FILTERING AND RECALIBRATION
    #########################################################
    out_vcfsnpr_name = group_name + ".cohort.snp.vcf"
    out_vcfindelr_name = group_name + ".cohort.indel.vcf"
    output_vcfsnpr_file = os.path.join(out_vcfr_dir, out_vcfsnpr_name)
    output_vcfindelr_file = os.path.join(out_vcfr_dir, out_vcfindelr_name)

    if os.path.isfile(output_vcfsnpr_file) and os.path.isfile(output_vcfindelr_file):
        logger.info(YELLOW + DIM + output_vcfsnpr_file + " EXIST\nOmmiting Variant Selection (Recall-Group) for group " + group_name + END_FORMATTING)
    else:
        logger.info(GREEN + "Selecting Variants (Recall-Group) in group " + group_name + END_FORMATTING)
        select_variants(output_vcfr_file, select_type='SNP') #select_variants(output_vcfr_file, select_type='INDEL')
        select_variants(output_vcfr_file, select_type='INDEL')

    #HARD FILTER VARIANTS 1/2 FOR RECALIBRATION #############
    #########################################################
    out_vcfhfsnpr_name = group_name + ".cohort.snp.hf.vcf"
    out_vcfhfindelr_name = group_name + ".cohort.indel.hf.vcf"
    output_vcfhfsnpr_file = os.path.join(out_vcfr_dir, out_vcfhfsnpr_name)
    output_vcfhfindelr_file = os.path.join(out_vcfr_dir, out_vcfhfindelr_name)


    if os.path.isfile(output_vcfhfsnpr_file) and os.path.isfile(output_vcfhfindelr_file):
        logger.info(YELLOW + DIM + output_vcfhfsnpr_file + " EXIST\nOmmiting Hard Filtering (Recall-Group) for group " + group_name + END_FORMATTING)
    else:
        logger.info(GREEN + "Hard Filtering Variants (Recall-Group) in group " + group_name + END_FORMATTING)
        hard_filter(output_vcfsnpr_file, select_type='SNP')
        hard_filter(output_vcfindelr_file, select_type='INDEL')

    #PASS FILTER VARIANTS 1/2 FOR RECALIBRATION #############
    #########################################################
    out_vcfhfsnppass_name = group_name + ".cohort.snp.hf.pass.vcf"
    #out_vcfhfindelpass_name = group_name + ".cohort.indel.hf.pass.vcf"
    output_vcfhfsnppass_file = os.path.join(out_vcfr_dir, out_vcfhfsnppass_name)
    #output_vcfhfindelpass_file = os.path.join(out_vcfr_dir, out_vcfhfindelpass_name)


    if os.path.isfile(output_vcfhfsnpr_file) and os.path.isfile(output_vcfhfsnppass_file):
        logger.info(YELLOW + DIM + output_vcfhfsnppass_file + " EXIST\nOmmiting PASS Filtering (Recall-Group) for group " + group_name + END_FORMATTING)
    else:
        logger.info(GREEN + "PASS Filtering Variants (Recall-Group) in group " + group_name + END_FORMATTING)
        select_pass_variants(output_vcfhfsnpr_file, nocall_fr=0.2)
        select_pass_variants(output_vcfhfindelr_file, nocall_fr=0.2)


        ######################################################################
        ##############START RECALIBRATION AND FINAL CALL######################
        ######################################################################

    logger.info("\n\n" + BLUE + BOLD + "STARTING RECALIBATION IN GROUP: " + group_name + END_FORMATTING)

    for r1_file, r2_file in zip(r1, r2):
        sample = extract_sample(r1_file, r2_file)

        args.sample = sample
        args.output = os.path.abspath(args.output)

        if sample in sample_list_F:

            sample_number = str(sample_list_F.index(sample) + 1)
            sample_total = str(len(sample_list_F))

            logger.info("\n" + WHITE_BG + "RECALIBRATION AND CALL ON SAMPLE: " + sample + " (" + sample_number + "/" + sample_total + ")" + END_FORMATTING)

            ##############START BAM RECALIBRATION############
            #################################################

            ################BQSR AND APPLY BQSR##################
            #####################################################
            out_bqsr_name = sample + ".bqsr.bam"
            output_bqsr_file = os.path.join(out_map_dir, out_bqsr_name)

            if os.path.isfile(output_bqsr_file):
                logger.info(YELLOW + DIM + output_bqsr_file + " EXIST\nOmmiting Recalibration for sample " + sample + END_FORMATTING)
            else:
                logger.info(GREEN + "Recalibration in sample " + sample + END_FORMATTING)
                recalibrate_bam(args)

            #HAPLOTYPE CALL 1/2 FOR HARD FILTERING AND RECALIBRATION
            #######################################################
            out_gvcf_name = sample + ".g.vcf"
            output_gvcf_file = os.path.join(out_gvcf_dir, out_gvcf_name)

            #args.input_bam = output_bqsr_file

            if os.path.isfile(output_gvcf_file):
                logger.info(YELLOW + DIM + output_gvcf_file + " EXIST\nOmmiting Haplotype Call for sample " + sample + END_FORMATTING)
            else:
                logger.info(GREEN + "Haplotype Calling in sample " + sample + END_FORMATTING)
                haplotype_caller(args, recalibrate=False, ploidy=2, bamout=False, forceactive=False)

    #ONCE ALL GVCF VARIANTS ARE CALLED, THEY ARE GATHERED AND FILTERED 
    # FOR FINAL CALLING
    ######################################################################
    ##############START GROUP CALLING FOR FINAL CALL######################
    ######################################################################
    group_name = args.output.split("/")[-1]
    logger.info("\n\n" + BLUE + BOLD + "STARTING JOINT CALL FOR FINAL CALL IN GROUP: " + group_name + END_FORMATTING + "\n")

    #CALL VARIANTS 2/2 FOR HARD FILTERING AND RECALIBRATION
    #######################################################
    out_gvcf_name = group_name + ".cohort.g.vcf"
    output_gvcf_file = os.path.join(out_gvcf_dir, out_gvcf_name)

    if os.path.isfile(output_gvcf_file):
        logger.info(YELLOW + DIM + output_gvcfr_file + " EXIST\nOmmiting GVCF Combination for group " + group_name + END_FORMATTING)
    else:
        logger.info(GREEN + "GVCF Combination in group " + group_name + END_FORMATTING)
        combine_gvcf(args, recalibrate=False, all_gvcf=args.enrich_gvcf)

    #CALL VARIANTS 2/2 FOR HARD FILTERING AND RECALIBRATION
    #######################################################
    out_vcf_name = group_name + ".cohort.raw.vcf"
    output_vcf_file = os.path.join(out_vcf_dir, out_vcf_name)

    if os.path.isfile(output_vcf_file):
        logger.info(YELLOW + DIM + output_vcf_file + " EXIST\nOmmiting Variant Calling (Group) for group " + group_name + END_FORMATTING)
    else:
        logger.info(GREEN + "Variant Calling (Group) in group " + group_name + END_FORMATTING)
        call_variants(args, recalibrate=False, group=True)

    #SELECT VARIANTS 2/2 FOR HARD FILTERING AND RECALIBRATION
    #########################################################
    out_vcfsnp_name = group_name + ".cohort.snp.vcf"
    out_vcfindel_name = group_name + ".cohort.indel.vcf"
    output_vcfsnp_file = os.path.join(out_vcf_dir, out_vcfsnp_name)
    output_vcfindel_file = os.path.join(out_vcf_dir, out_vcfindel_name)

    if os.path.isfile(output_vcfsnp_file) and os.path.isfile(output_vcfindel_file):
        logger.info(YELLOW + DIM + output_vcfsnp_file + " EXIST\nOmmiting Variant Selection (Group) for group " + group_name + END_FORMATTING)
    else:
        logger.info(GREEN + "Selecting Variants (Group) in group " + group_name + END_FORMATTING)
        select_variants(output_vcf_file, select_type='SNP')
        select_variants(output_vcf_file, select_type='INDEL')

    #HARD FILTER VARIANTS 2/2 FOR RECALIBRATION #############
    #########################################################
    out_vcfhfsnp_name = group_name + ".cohort.snp.hf.vcf"
    out_vcfhfindel_name = group_name + ".cohort.indel.hf.vcf"
    output_vcfhfsnp_file = os.path.join(out_vcf_dir, out_vcfhfsnp_name)
    output_vcfhfindel_file = os.path.join(out_vcf_dir, out_vcfhfindel_name)


    if os.path.isfile(output_vcfhfsnp_file) and os.path.isfile(output_vcfhfindel_file):
        logger.info(YELLOW + DIM + output_vcfhfsnp_file + " EXIST\nOmmiting Hard Filtering (Group) for group " + group_name + END_FORMATTING)
    else:
        logger.info(GREEN + "Hard Filtering Variants (Group) in group " + group_name + END_FORMATTING)
        hard_filter(output_vcfsnp_file, select_type='SNP')
        hard_filter(output_vcfindel_file, select_type='INDEL')


    #PASS FILTER VARIANTS 2/2 FOR RECALIBRATION #############
    #########################################################
    out_vcfhfcombined_name = group_name + ".cohort.combined.hf.vcf"
    output_vcfhfcombined_file = os.path.join(out_vcf_dir, out_vcfhfcombined_name)


    if os.path.isfile(output_vcfhfcombined_file):
        logger.info(YELLOW + DIM + output_vcfhfcombined_file + " EXIST\nOmmiting combination for group " + group_name + END_FORMATTING)
    else:
        logger.info(GREEN + "Combining both vcf SNP and INDEL in group " + group_name + END_FORMATTING)
        combine_vcf(output_vcfhfsnp_file, output_vcfhfindel_file, name_out=False)

    if args.all_cohort == True:
        split_vcf_saples(output_vcfhfcombined_file, sample_list=False, nocall_fr=args.maxnocallfr)
    else:
        split_vcf_saples(output_vcfhfcombined_file, sample_list=sample_list_F, nocall_fr=args.maxnocallfr)
    ###########################################################################
    ###########################################################################
    ###########################################################################

    logger.info(GREEN + "Determinind highly heterozygous and poorly genotyped regions in " + group_name + END_FORMATTING)
    highly_hetz_to_bed(output_vcfhfcombined_file, "highly_hetz", reference="CHROM", nocall_fr=0.5)
    non_genotyped_to_bed(output_vcfhfcombined_file, "non_genotyped", reference="CHROM", nocall_fr=0.5)

    for r1_file, r2_file in zip(r1, r2):
        sample = extract_sample(r1_file, r2_file)
        args.sample = sample
        args.output = os.path.abspath(args.output)

        if sample in sample_list_F:

            logger.info("\n" + WHITE_BG + "FINAL FILTERING IN SAMPLE " + sample + END_FORMATTING)

            ################FINAL VCF FILTERING##################
            #####################################################
            out_final_name = sample + ".combined.hf.ALL.final.vcf"
            in_final_name = sample + ".combined.hf.vcf"
            output_final_vcf = os.path.join(out_vcf_dir, out_final_name)
            in_final_vcf = os.path.join(out_vcf_dir, in_final_name)

            if os.path.isfile(output_final_vcf):
                logger.info(YELLOW + DIM + output_final_vcf + " EXIST\nOmmiting Final filter for sample " + sample + END_FORMATTING)
            else:
                logger.info(GREEN + "Final filter in sample " + sample + END_FORMATTING)
                vcf_consensus_filter(in_final_vcf, distance=1, AF=0.75, QD=15, window_10=3, dp_limit=8, dp_AF=10, AF_dp=0.80,
                highly_hetz=highly_hetz_bed, 
                non_genotyped=non_genotyped_bed, 
                poorly_covered=poorly_covered_bed, 
                var_type="SNP")

    #DETEMINING MIXED ORIGIN IN GROUP######################
    #######################################################
    output_vcfstat_file = os.path.join(out_table_dir, "vcf_stat.tab")
    if os.path.isfile(output_vcfstat_file):
        logger.info("\n" + YELLOW + DIM + output_vcfstat_file + " EXIST\nOmmiting Mixed search in group " + group_name + END_FORMATTING)
        samples_mixed = []
    else:
        logger.info(GREEN + "Finding Mixed samples in " + group_name + END_FORMATTING)
        samples_mixed = vcf_stats(out_table_dir, distance=15, quality=10)

    if len(samples_mixed) > 0:
        logger.info("\n" + YELLOW + BOLD + "There are mixed sample(s): " + "\n"\
            + ",".join(samples_mixed) + END_FORMATTING + "\n")
        remove_low_covered_mixed(args.output, samples_mixed, "Mixed")
        #Remove sample from the list of filtered samples
        ################################################
        for samples_to_remove in samples_mixed:
            sample_list_F.remove(samples_to_remove)
    else:
        logger.info("\n" + YELLOW + BOLD + "No mixed samples have been detected" + "\n")

    logger.info("\n\n" + MAGENTA + BOLD + "VARIANT CALL FINISHED IN GROUP: " + group_name + END_FORMATTING + "\n")

    #######################################################################################################################################
    #################################END OF VARIANT CALLING################################################################################
    #######################################################################################################################################
    tuberculosis = False
    if tuberculosis == True:
        logger.info("\n\n" + BLUE + BOLD + "STARTING ANNOTATION IN GROUP: " + group_name + END_FORMATTING + "\n")

        for root, _, files in os.walk(out_vcf_dir):
            for name in files:
                filename = os.path.join(root, name)
                output_path = os.path.join(out_annot_dir, name)
                if filename.endswith("combined.hf.vcf"):
                    sample = name.split(".")[0]
                    if sample in sample_list_F:
                        #ANNOTATION -AUTO AND SPECIFIC- ###################
                        ###################################################
                        out_annot_name = sample + ".combined.hf.annot.tsv"
                        output_annot_file = os.path.join(out_annot_dir, out_annot_name)

                        if os.path.isfile(output_annot_file):
                            logger.info(YELLOW + DIM + output_annot_file + " EXIST\nOmmiting Annotation for sample " + sample + END_FORMATTING)
                        else:
                            logger.info(GREEN + "Annotating snps in sample " + sample + END_FORMATTING)
                            replace_reference(filename, output_path)
                            snpeff_annotation(args, output_path, database=args.snpeff_database)
                            #Handle output vcf file from SnpEff annotation
                            vcf_path = (".").join(output_path.split(".")[:-1])
                            annot_vcf = vcf_path + ".annot"
                            #This function add SPECIFIC anotation
                            if args.annot_bed:
                                final_annotation(annot_vcf, *args.annot_bed)
                            else:
                                final_annotation(annot_vcf)


        logger.info("\n\n" + MAGENTA + BOLD + "ANNOTATION FINISHED IN GROUP: " + group_name + END_FORMATTING + "\n")
    else:
        logger.info("NO TB Selected, snpEff won't be executed")




    logger.info("\n\n" + BLUE + BOLD + "STARTING COMPARISON IN GROUP: " + group_name + END_FORMATTING + "\n")

    check_create_dir(out_compare_dir)
    folder_compare = today + "_" + group_name
    path_compare = os.path.join(out_compare_dir, folder_compare)
    check_create_dir(path_compare)
    full_path_compare = os.path.join(path_compare, group_name)

    

    #ddtb_add(out_vcf_dir, full_path_compare)
    ddtb_add(out_vcf_dir, full_path_compare, recalibrate=args.output)

    compare_snp_matrix = full_path_compare + ".revised.tsv"
    
    ddtb_compare(compare_snp_matrix)

    logger.info("\n\n" + MAGENTA + BOLD + "COMPARING FINISHED IN GROUP: " + group_name + END_FORMATTING + "\n")


    if args.noclean == True:
        logger.info("\n\n" + BLUE + BOLD + "STARTING CLEANING IN GROUP: " + group_name + END_FORMATTING + "\n")
        clean_unwanted_files(args)
    else:
        logger.info("No cleaning was requested")

    
    """
    logger.info("\n\n" + MAGENTA + BOLD + "#####END OF PIPELINE SNPTB#####" + END_FORMATTING + "\n")

if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        logger.exception(e)
        raise
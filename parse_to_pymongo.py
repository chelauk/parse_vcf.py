import sys                         # for moving files around 
import re                          # to create regular expressions 
from pymongo import MongoClient    # for adding to mongodb
import pysam                       # for parsing a tabix indexed vcf file
from collections import defaultdict

'''
I may or may not include flask in in this parser
import json
from flask import Flask, request, session, g, redirect, url_for, abort, render_template, flash, jsonify, send_from_directory
'''

client = MongoClient()           # create a link to mongodb
db = client.vcf_db               # an a variable to represent the monogodb database, either existing or creates it
variants=db.variant_collection   # a variable to represent the collection either existing or creates it
# for testing clear variants everytime
variants.remove({ })             # while tesing
# a function to parse the vcf-file header to get samples
vcf = pysam.Tabixfile("/Volumes/DATA/vcfs/gosgene_ann_vcf_db.vcf.gz") # an variable to hold the vcf file
vcf_header = vcf.header
body = vcf.fetch(1, 0, 103000)

# create annotions regular expression

ann_pattern = re.compile('ANN=*')
allelle_count_pattern = re.compile('AC=')

def parse_header_for_samples():
    temp_list = list( vcf.header )[-1] # last item of iteritem
    temp_list = temp_list.split('\t')[9:]
    return temp_list

def parse_body_line():
    for line in body:
        count = 0
        body_list = line.split('\t')
        sample_annotations=body_list[7]
        sample_annotations=sample_annotations.split(';')  # variant effect predictor annotaions (taking the first set)
        for item in sample_annotations:
            if allelle_count_pattern.match(item):
                allelle_count = item
                allelle_count = allelle_count.strip('AC=')
        for item in sample_annotations:
            if ann_pattern.match(item):
                sample_annotations = item
        sample_annotations=sample_annotations.split('|')
        sample_genotypes = body_list[8]
        genotype_fields = sample_genotypes.split(':')
        chr_pos_ref_alt = str(body_list[0]) + '_' + str(body_list[1]) + '_' + str(body_list[3]) + '_' + str(body_list[4])
        position_dict = defaultdict(list)
        position_dict[chr_pos_ref_alt] = [ body_list[0] , body_list[1] , body_list[3], body_list[4] ]
        sample_dict = defaultdict(list)
        sample_info = body_list[9:]
        sample_gt_info = []   #       selecting those samples with a genotype at this position
        for i in sample_info:
            affected_sample=0
            count += 1
            first_character = i[0:1]
            if first_character is not '.':
                affected_sample += 1
                sample_genotype_breakdown = i.split(':')
                #print genotype_fields
                if 'DP' not in genotype_fields:
                    genotype_fields.append('DP')
                    sample_genotype_breakdown.append(1)
                if 'AD' not in genotype_fields:
                    genotype_fields.append('AD')
                    sample_genotype_breakdown.append(1)
                
                genotype_dict = dict(zip(genotype_fields,sample_genotype_breakdown))
                genotype_list = {}
                sample_genotype_breakdown = [genotype_dict['GT'], genotype_dict['AD']]
                this_sample = sample_list[count].replace('.', '_')#print genotype_dict['GT'] + '\t' + str(genotype_dict['DP']) + '\t' + str(genotype_dict['AD']) + "\t"  + sample_list[count]
                genotype_count = 0 # variable for generating average depth
                cumulative_depth = 0 # variable for generating average depth
        
                if chr_pos_ref_alt == str(body_list[0]) + '_' + str(body_list[1]) + '_' + str(body_list[3]) + '_' + str(body_list[4]):
                    try:
                        position_dict[chr_pos_ref_alt][4]
                    except IndexError:
                        position_dict[chr_pos_ref_alt].append(sample_annotations[1])
                    try:
                        position_dict[chr_pos_ref_alt][5]
                    except IndexError:
                        position_dict[chr_pos_ref_alt].append(sample_annotations[3])                                        
                    try:    
                        position_dict[chr_pos_ref_alt][6]
                    except IndexError:
                        position_dict[chr_pos_ref_alt].append(sample_annotations[10])
                    try:    
                        position_dict[chr_pos_ref_alt][7]
                    except IndexError:
                        position_dict[chr_pos_ref_alt].append(allelle_count)
                    
                    
                    
                    
                    try: 
                        sample_dict[chr_pos_ref_alt][0]
                    except IndexError:
                        sample_dict[chr_pos_ref_alt].append(genotype_list)
                    sample_dict[chr_pos_ref_alt][0][this_sample] = sample_genotype_breakdown
                    try:
                        sample_depths = sample_genotype_breakdown[1].split(',')
                        sample_depths = map(int, sample_depths)
                    except AttributeError:
                        sample_depths = sample_genotype_breakdown[1]
                    try:
                        tot_sample_depth = sum(sample_depths)
                    except TypeError:
                        tot_sample_depths = sample_depths
                    cumulative_depth += tot_sample_depth
                    #print "depth " + cumulative_depth
                    genotype_count += 1
               
        variants.insert( {"_id":chr_pos_ref_alt,
                          "variant" : position_dict[chr_pos_ref_alt],
                          "sample" : sample_dict[chr_pos_ref_alt]
                          })
                #print str(body_list[0]) + '\t' + str(body_list[1]) + '\t' + str(body_list[2]) + '\t' + str(body_list[3]) + '\t' \
                #+ str(body_list[4]) + '\t' + i + "\t" + genotype_dict['GT'] + '\t' + str(genotype_dict['DP']) + '\t' + str(genotype_dict['AD']) + "\t"  + sample_list[count]
    

sample_list = parse_header_for_samples()

parse_body_line()

for i in variants.find({'variant': { "$all":['1','17594']}}): # to use mongo $all it needs to be quotes
    print i
#!/usr/bin/env python
# coding: utf-8

import re
import os
import sys
import time
import logging
from Bio import SeqIO, Entrez, Data
from .config import CommonNamesDict, MTColors
from Bio.SeqFeature import CompoundLocation, SimpleLocation

class Feature:
    def __init__(self, name, location, type, color, join=None, mtgenome=None, accession=None, file=None, topology=None, partition=None):
        self.name = name
        self.location = location
        self.type = type
        self.color = color
        self.join = join
        self.mtgenome = mtgenome
        self.accession = accession
        self.file = file
        self.topology = topology
        self.partition = partition
        
    def __repr__(self):
            return str(self.name)
    def __str__(self):
            return self.__repr__()

def search_name(gene_name):
    gene_name = gene_name.upper()
    status = False
    if CommonNamesDict.get(gene_name) == None:
        for p in CommonNamesDict:
            if bool(re.search("^"+re.escape(p)+".*", gene_name)):
                status = True
                return CommonNamesDict.get(p)
    if status == False:
        return gene_name
    else:
        return CommonNamesDict.get(gene_name)


def is_repeat(features, location, name=None):
    status = False
    for i in features:
        if i.location == location:
            status = True

        if i.type=="CDS" and i.name == name:
            if location.strand == i.location.strand:
                if (location.start <= i.location.start and location.end >= i.location.end) or \
                (i.location.start <= location.start and i.location.end >= location.end):
                    status = True
    return status

def rotate_seq(seq, index):
    k = len(seq) - index
    k = k % len(seq)
    return seq[-k:] + seq[:-k]

def get_genbank_from_ncbi(accession):
    Entrez.email = "thecgs001.foxmail.com"
    gb_text = Entrez.efetch(db="Nucleotide", id=accession, rettype='gb')
    return gb_text
    
def get_species_name(string, abbr=True):
    species = re.sub("_", " ", string)
    if abbr:
        if bool(re.search(' x ', species)):
            tmp = species.split(' x ')
            species1 = tmp[0].split(' ')[0][0] + '. ' +  tmp[0].split(' ')[1]
            species2 = tmp[1].split(' ')[0][0] + '. ' +  tmp[1].split(' ')[1]
            species = species1 + ' x ' + species2
        else:
            if len(species.split(' ')) == 2:
                species = species.split(' ')[0][0] + '. ' + species.split(' ')[1]
            if len(species.split(' ')) == 3:
                species = species.split(' ')[0][0] + '. ' + species.split(' ')[1] + ' ' + ' '.join(species.split(' ')[2:])
    return species

def reinit_features(features, start = "tRNA-Phe", force_reoriented=False):
    logger = logging.getLogger(__name__) 
    logger.setLevel(logging.DEBUG)                  
    features_new = [features[0]]
    
    value = None
    if isinstance(start, str):
        for feature in features[1:]:
            if feature.name == start:
                value = feature.location.start
                break
                
    if value == None:
        logger.warning(features[0].file+" does not contain" + start + ". Therefore, the rotation operation was not performed. Please check.")
        return features
        
    if features[0].topology == "linear" and force_reoriented == False:
        logger.warning(features[0].file+" is linear topology. Reoriention is disabled by default, so no rotation operation was performed. To force reoriention, use the force_reoriented=True parameter.")
        return features
        
    _status = True
    for index, feature in enumerate(features[1:]):
        #print(feature, feature.location, feature.name)
        if feature.name == start:
            _status = False
        if _status:
            location = SimpleLocation(start=feature.location.start + len(features[0].location)-value,
                                      end=feature.location.end + len(features[0].location)-value,
                                      strand=feature.location.strand)
            features_new.append(Feature(feature.name, location, feature.type, feature.color, feature.join))
        else:
            location = SimpleLocation(start=feature.location.start-value,
                                      end=feature.location.end-value,
                                      strand=feature.location.strand)
            features_new.append(Feature(feature.name, location, feature.type, feature.color, feature.join))
            
    tmp = [features_new[0]]
    tmp.extend(sorted(features_new[1:], key=lambda x:x.location.start))
    f = Feature(features[0].name, features[0].location, features[0].type, features[0].color, features[0].join, rotate_seq(features[0].mtgenome, value), accession=features[0].accession, file=features[0].file, topology=features[0].topology, partition=features[0].partition)
    features_new = [f]
    index = 0
    for feature in tmp[1:]:
        if (features_new[index].name != feature.name) or (features_new[index].join==None) or (feature.join==None):
            features_new.append(feature)
            index += 1
        else:
            location = SimpleLocation(start=features_new[index].location.start,
                                      end=feature.location.end,
                                      strand=feature.location.strand)
            features_new[index].location = location
            features_new[index].join = None
            
    return features_new


def get_type(genename):
    if genename in ['ND1', 'ND2', 'ND3', 'ND4L', 'ND4', 'ND5', 'ND6', 'COX1', 'COX2', 'COX3', 'ATPase6', 'ATPase8', 'Cytb']:
        return "CDS"
    elif 'rRNA' in genename:
        return "rRNA"
    elif 'tRNA' in genename:
        return "tRNA"
    else:
        return "CDS"
        
def get_features(file, abbr=False, colors=None, isfilename2species=False, start=None, force_reoriented=False):
    """
    Descripton:
        Parses a genbank file and returns a list of feature classes.
    
    Parameters：
        file: {str} one genbankfile or NCBI accession ID.
        abbr: {bool} whether to abbreviate species names.
        isfilename2species: {bool} whether filename convert to species.
        start: {None, str} initial feature, such as, ND1, ND2, ND3, ND4, ND4L, ND5, ND6,
                     COX1, COX2, COX3, ATPase6, ATPase8, Cytb, tRNA-His, tRNA-Pro,
                     tRNA-Thr, tRNA-Trp, tRNA-Met, tRNA-Asp, tRNA-Ala, tRNA-Gln,
                     tRNA-Ile, tRNA-Arg, tRNA-Tyr, tRNA-Phe, tRNA-Lys, tRNA-Gly,
                     tRNA-Asn, tRNA-Leu, tRNA-Glu, tRNA-Val, tRNA-Cys, tRNA-Ser,
                     12S rRNA, 16S rRNA, D-loop. default=None.   
        colors: {str, dict} themes such as, Chen, Tan, ogdraw, mitofish,
                            mitofish1, mitoz,  gggenes, chloroplot, grey, igv.
        force_reoriented: {bool} force-reoriendted linear mtgenome.
    """
    if colors == None:
        colors = MTColors['MITOFISH']
    elif isinstance(colors, str):
        colors = MTColors.get(colors.upper(), 'MITOFISH')
    elif isinstance(colors, dict):
        pass
            
    features = []
    unknown_locations = {}
    
    if os.path.exists(file):
        handle = SeqIO.parse(file, 'genbank')
    else:
        handle = SeqIO.parse(get_genbank_from_ncbi(file), 'genbank')
        
    #print(handle)
    for record in handle:
        if 'data_file_division' in record.annotations:
            partition = record.annotations['data_file_division']
        else:
            partition = "UNA"
            
        topology = record.annotations['topology'].lower()
        mtgenome = record.seq.upper()
        accession = record.id
        
        for i in record.features:
            if len(list(i.qualifiers.values())) != 0:
                if "product" in i.qualifiers:
                    gene_name = i.qualifiers['product'][0]
                    gene_name =  search_name(gene_name)
                    
                    #print("gene_product", gene_name)
                    if gene_name.upper() not in CommonNamesDict and "gene" in i.qualifiers:
                        gene_name = i.qualifiers['gene'][0]
                        gene_name =  search_name(gene_name)
                        
                        if gene_name.upper() not in CommonNamesDict and "note" in i.qualifiers:
                            gene_name = i.qualifiers['note'][0]
                            gene_name =  search_name(gene_name)
                            if gene_name.upper() not in CommonNamesDict and i.type.upper() in CommonNamesDict:
                                gene_name =  search_name(i.type)
                        
                elif "gene" in i.qualifiers:
                    gene_name = i.qualifiers['gene'][0]
                    gene_name =  search_name(gene_name)
                    #print(i, gene_name)
                    if gene_name.upper() not in CommonNamesDict and "note" in i.qualifiers:
                        #print(i, gene_name)
                        gene_name = i.qualifiers['note'][0]
                        gene_name =  search_name(gene_name)
                        if gene_name.upper() not in CommonNamesDict and i.type.upper() in CommonNamesDict:
                            gene_name =  search_name(i.type)
                
                elif "note" in i.qualifiers:
                    #print(i.qualifiers)
                    gene_name = i.qualifiers["note"][0]
                    gene_name =  search_name(gene_name)
                    #print(i, gene_name)
                    if gene_name.upper() not in CommonNamesDict and i.type.upper() in CommonNamesDict:
                        gene_name =  search_name(i.type)
                                
                elif "organism" in i.qualifiers:
                    gene_name = i.qualifiers["organism"][0]
                    #gene_name =  search_name(gene_name)
                    
                else:
                    #print(i.qualifiers)
                     continue
            else:
                gene_name = i.type
                gene_name =  search_name(gene_name)
                
            gene_name = CommonNamesDict.get(gene_name.upper(), gene_name)
            ## repalce
            if gene_name.upper() not in CommonNamesDict:
                if str(i.location) not in unknown_locations:
                    unknown_locations[str(i.location)] = (gene_name, False, i.location)
            else:
                if str(i.location) in unknown_locations:
                    unknown_locations[str(i.location)] = (gene_name, True, i.location)
                    
            if i.type ==  "source":
                if isfilename2species:
                    species_name = os.path.splitext(os.path.basename(file))[0]
                else:
                    species_name = i.qualifiers["organism"][0]
                    
                species_name = get_species_name(species_name, abbr=abbr)
                features.append(Feature(name=species_name, location=i.location, type=i.type, color=colors.get('source', colors.get('Other genes', 'gray')),
                                        mtgenome=mtgenome, accession=accession, file=file, topology=topology, partition=partition))

            elif i.type in ['rRNA', 'tRNA', 'D_loop', 'D-loop']:
                if gene_name in ['tRNA-His', 'tRNA-Pro', 'tRNA-Thr', 'tRNA-Trp', 'tRNA-Met', 'tRNA-Asp', 'tRNA-Ala', 'tRNA-Gln',
                                 'tRNA-Ile', 'tRNA-Arg', 'tRNA-Tyr', 'tRNA-Phe', 'tRNA-Lys', 'tRNA-Gly', 'tRNA-Asn', 'tRNA-Leu',
                                 'tRNA-Glu', 'tRNA-Val', 'tRNA-Cys', 'tRNA-Ser', '12S rRNA', '16S rRNA', "D-loop"]:
                    if isinstance(i.location, CompoundLocation):
                        for location in i.location.parts:
                            if not is_repeat(features=features, location=location):
                                features.append(Feature(name=gene_name, location=location, type=i.type, color=colors.get(gene_name, colors.get('Other genes', 'gray')),join=i.location))
                    else:
                        if not is_repeat(features=features, location=i.location):
                            features.append(Feature(name=gene_name, location=i.location, type=i.type, color=colors.get(gene_name, colors.get('Other genes', 'gray'))))                   
                    
            elif i.type in ['CDS', 'gene']:
                if gene_name in ['ND1', 'ND2', 'ND3', 'ND4L', 'ND4', 'ND5', 'ND6', 'COX1', 'COX2', 'COX3', 'ATPase6', 'ATPase8', 'Cytb']:
                    if isinstance(i.location, CompoundLocation):
                        for location in i.location.parts:
                            if not is_repeat(features=features, location=location, name=gene_name):
                                features.append(Feature(name=gene_name, location=location, type="CDS", color=colors.get(gene_name, colors.get('Other genes', 'gray')),join=i.location))
                    else:
                        if not is_repeat(features=features, location=i.location, name=gene_name):
                            features.append(Feature(name=gene_name, location=i.location, type="CDS", color=colors.get(gene_name, colors.get('Other genes', 'gray'))))
                
                elif 'tRNA' in gene_name:
                    if isinstance(i.location, CompoundLocation):
                        for location in i.location.parts:
                            if not is_repeat(features=features, location=location, name=gene_name):
                                features.append(Feature(name=gene_name, location=location, type="tRNA", color=colors.get(gene_name, colors.get('Other genes', 'gray')),join=i.location))
                    else:
                        if not is_repeat(features=features, location=i.location, name=gene_name):
                            features.append(Feature(name=gene_name, location=i.location, type="tRNA", color=colors.get(gene_name, colors.get('Other genes', 'gray'))))
                            
                elif gene_name in ['12S rRNA', '16S rRNA']:
                    if isinstance(i.location, CompoundLocation):
                        for location in i.location.parts:
                            if not is_repeat(features=features, location=location, name=gene_name):
                                features.append(Feature(name=gene_name, location=location, type="rRNA", color=colors.get(gene_name, colors.get('Other genes', 'gray')),join=i.location))
                    else:
                        if not is_repeat(features=features, location=i.location, name=gene_name):
                            features.append(Feature(name=gene_name, location=i.location, type="rRNA", color=colors.get(gene_name, colors.get('Other genes', 'gray'))))
                else: #ORF
                    if isinstance(i.location, CompoundLocation):
                        for location in i.location.parts:
                            if not is_repeat(features=features, location=location, name=gene_name):
                                features.append(Feature(name=gene_name, location=location, type="CDS", color=colors.get(gene_name, colors.get('Other genes', 'gray')),join=i.location))
                    else:
                        if not is_repeat(features=features, location=i.location, name=gene_name):
                            features.append(Feature(name=gene_name, location=i.location, type="CDS", color=colors.get(gene_name, colors.get('Other genes', 'gray'))))

                            
            elif i.type in ['misc_feature', 'repeat_region']:
                if gene_name in ['tRNA-His', 'tRNA-Pro', 'tRNA-Thr', 'tRNA-Trp', 'tRNA-Met', 'tRNA-Asp', 'tRNA-Ala', 'tRNA-Gln',
                                 'tRNA-Ile', 'tRNA-Arg', 'tRNA-Tyr', 'tRNA-Phe', 'tRNA-Lys', 'tRNA-Gly', 'tRNA-Asn', 'tRNA-Leu',
                                 'tRNA-Glu', 'tRNA-Val', 'tRNA-Cys', 'tRNA-Ser']:
                    if isinstance(i.location, CompoundLocation):
                        for location in i.location.parts:
                            if not is_repeat(features=features, location=location, name=gene_name):
                                features.append(Feature(name=gene_name, location=location, type="tRNA", color=colors.get(gene_name, colors.get('Other genes', 'gray')),join=i.location))
                    else:
                        if not is_repeat(features=features, location=i.location, name=gene_name):
                            features.append(Feature(name=gene_name, location=i.location, type="tRNA", color=colors.get(gene_name, colors.get('Other genes', 'gray'))))
                            
                elif gene_name in ['12S rRNA', '16S rRNA']: #, "D-loop"
                    if isinstance(i.location, CompoundLocation):
                        for location in i.location.parts:
                            if not is_repeat(features=features, location=location, name=gene_name):
                                features.append(Feature(name=gene_name, location=location, type="rRNA", color=colors.get(gene_name, colors.get('Other genes', 'gray')),join=i.location))
                    else:
                        if not is_repeat(features=features, location=i.location, name=gene_name):
                            features.append(Feature(name=gene_name, location=i.location, type="rRNA", color=colors.get(gene_name, colors.get('Other genes', 'gray'))))
                        
                elif gene_name in ["D-loop"]:
                    if isinstance(i.location, CompoundLocation):
                        for location in i.location.parts:
                            if not is_repeat(features=features, location=location, name=gene_name):
                                features.append(Feature(name=gene_name, location=location, type="D-loop", color=colors.get(gene_name, colors.get('Other genes', 'gray')),join=i.location))
                    else:
                        if not is_repeat(features=features, location=i.location, name=gene_name):
                            features.append(Feature(name=gene_name, location=i.location, type="D-loop", color=colors.get(gene_name, colors.get('Other genes', 'gray'))))
                        
                        
                elif gene_name in ['ND1', 'ND2', 'ND3', 'ND4L', 'ND4', 'ND5', 'ND6', 'COX1', 'COX2', 'COX3', 'ATPase6', 'ATPase8', 'Cytb']:
                    if isinstance(i.location, CompoundLocation):
                        for location in i.location.parts:
                            if not is_repeat(features=features, location=location, name=gene_name):
                                features.append(Feature(name=gene_name, location=location, type="CDS", color=colors.get(gene_name, colors.get('Other genes', 'gray')),join=i.location))
                    else:
                        if not is_repeat(features=features, location=i.location, name=gene_name):
                            features.append(Feature(name=gene_name, location=i.location, type="CDS", color=colors.get(gene_name, colors.get('Other genes', 'gray'))))
            else:
                pass


    join_locations = []
    features_tmp = [features[0]]
    for i in features[1:]:
        if i.join !=None:
            if str(i.join) in unknown_locations:
                if unknown_locations[str(i.join)][1] == True:
                    if str(i.join) not in join_locations:
                        join_locations.append(str(i.join))
                        for location in i.join.parts:
                            features_tmp.append(Feature(name=unknown_locations[str(i.join)][0],
                                                        location=location, 
                                                        type=get_type(unknown_locations[str(i.join)][0]),
                                                        color=colors.get(unknown_locations[str(i.join)][0],
                                                                         colors.get('Other genes', 'gray')),join=i.join))
                    else:
                        pass
                else:
                    features_tmp.append(i)
            else:
                features_tmp.append(i)
        else:
            if str(i.location) in unknown_locations:
                if unknown_locations[str(i.location)][1] == True:                
                    features_tmp.append(Feature(name=unknown_locations[str(i.location)][0],
                                                location=i.location,
                                                type=get_type(unknown_locations[str(i.location)][0]),
                                                color=colors.get(unknown_locations[str(i.location)][0],
                                                                 colors.get('Other genes', 'gray'))))                           
                else:
                    features_tmp.append(i)
            else:
                features_tmp.append(i)

    
    #print(features)
    features = features_tmp
    if features == []:
        logger = logging.getLogger(__name__) 
        logger.setLevel(logging.DEBUG)
        logger.warning(file + " Gene ID error!")
        sys.exit()
    
    res = [features[0]]
    
    res.extend(sorted(features[1:], key=lambda x:x.location.start))
    
    if res[0].topology == "linear":
        res.append(Feature(name="Gap", type="Gap", color='black',
                           join=None, mtgenome=None, accession=None, file=None, topology=None,
                           location=SimpleLocation(start=len(res[0].mtgenome)-2, end=len(res[0].mtgenome)+1, strand=0))
                  )	              
    if start !=None:
        res = reinit_features(res, start = start, force_reoriented=force_reoriented)
    return res
    
def tidy_genbank(file, output=None, isfilename2species=False, start=None, table=2, force_reoriented=False, partition="inherit"):
    """
    Descripton:
        Use PyVAMR's powerful GenBank parser to reorganize the GenBank 
        and generate a new GenBank file.
    
    Parameters：
        file: {str} a genbankfile or NCBI accession ID.
        tabe: {int} codon tables. such as 1-6, 9-16, 21-33.
        start: {None, str} initial feature, such as, ND1, ND2, ND3, ND4, ND4L, ND5, ND6,
                     COX1, COX2, COX3, ATPase6, ATPase8, Cytb, tRNA-His, tRNA-Pro,
                     tRNA-Thr, tRNA-Trp, tRNA-Met, tRNA-Asp, tRNA-Ala, tRNA-Gln,
                     tRNA-Ile, tRNA-Arg, tRNA-Tyr, tRNA-Phe, tRNA-Lys, tRNA-Gly,
                     tRNA-Asn, tRNA-Leu, tRNA-Glu, tRNA-Val, tRNA-Cys, tRNA-Ser,
                     12S rRNA, 16S rRNA, D-loop. default=None.
        output: {str} a path of genbank output file.
        partition: {str} a data file division. inherit: inherit, PRI: primate, ROD: rodent, MAM: mammal, VRT: vertebrate, INV: invertebrate,
                    PLN: plant, BCT: bacterial, VRL: viral, PHG: bacteriophage, SYN: synthetic,
                    UNA: unannotated, ENV: environmental sample.
    """
    product = {'ND1': 'NADH dehydrogenase subunit 1',
               'ND2': 'NADH dehydrogenase subunit 2',
               'ND3': 'NADH dehydrogenase subunit 3',
               'ND4L':'NADH dehydrogenase subunit 4L',
               'ND4': 'NADH dehydrogenase subunit 4',
               'ND5': 'NADH dehydrogenase subunit 5',
               'ND6': 'NADH dehydrogenase subunit 6',
               'COX1':'cytochrome c oxidase subunit 1',
               'COX2':'cytochrome c oxidase subunit 2',
               'COX3':'cytochrome c oxidase subunit 3',
               'ATPase6':'ATP synthase F0 subunit 6',
               'ATPase8':'ATP synthase F0 subunit 8',
               'Cytb':'cytochrome b', 
               '12S rRNA':'12S ribosomal RNA',
               '16S rRNA':'16S ribosomal RNA'}
    
    def textwrap_fill(s, w=58, sep='\n                     '):
        r = ''
        for i ,j in zip(range(0, len(s), w), range(w, len(s), w)):
            r +=s[i:j]+sep
        if len(s)%w==0:
            r +=s[-w:]
        else:
             r +=s[len(s)%w*-1:]
        return r

    def get_translation_string(feature, mtgenome, table):
        CDS = feature.location.extract(mtgenome)
        pep = str(CDS.translate(table=table))
        if str(CDS[0:3]) in Data.CodonTable.unambiguous_dna_by_id[2].start_codons:
            pep = "M" + pep[1:]
        if pep[-1] == "*":
            pep = pep[:-1]
        return textwrap_fill('/translation="'+pep+'"', w=58)

    def format_mtgenome(seq, line_length=60, group_length=10):
        formatted_lines = []
        for i in range(0, len(seq), line_length):
            line_start = i + 1
            line = seq[i:i+line_length]
            groups = [line[j:j+group_length] for j in range(0, len(line), group_length)]
            grouped_line = ' '.join(groups)
            formatted_lines.append(f"{line_start:>9d} {grouped_line}")    
        return '\n'.join(formatted_lines)

    features =  get_features(file, isfilename2species=isfilename2species, start=start, force_reoriented=force_reoriented)
    
    features_tmp = []
    tmp = []
    for feature in features:
        if feature.join == None:
            features_tmp.append(feature)
        else:
            if feature.join not in tmp:
                features_tmp.append(feature)
                tmp.append(feature.join)
    
    features = features_tmp
    organism = features[0].name
    genome_len = len(features[0].mtgenome)
    
    #print(features)
    
    if partition == "inherit":
        partition = features[0].partition
        
    gb_text = f"""LOCUS       {organism.replace(' ', '_')}                {genome_len} bp    DNA     {features[0].topology}     ​​{partition} {time.strftime("%d-%b-%Y", time.localtime()).upper()}
DEFINITION  .
ACCESSION   .
VERSION     .
KEYWORDS    .
SOURCE      mitochondrion {organism}
  ORGANISM  {organism}
            Unclassified.
REFERENCE   1  (bases 1 to {genome_len})
  AUTHORS   Chen, G.
  TITLE     PyVAM: A Python package for visualizing animal mitochondrial.
  JOURNAL   Unpublished
  TITLE     Direct Submission
FEATURES             Location/Qualifiers
     source          1..{genome_len}
                     /organism="{organism}"
                     /organelle="mitochondrion"
                     /mol_type="genomic DNA"
"""

    for feature in features[1:]:
        if feature.location.strand == 1:
            if feature.join == None:
                pos = f"{str(feature.location.start+1)}..{str(feature.location.end)}"
            else:
                pos = "join(" + ','.join( f"{i.start+1}..{i.end}"  for i in features[1].join.parts) + ")"
        elif feature.location.strand == -1:
            if feature.join == None:
                #print(feature.join)
                pos = f"complement({str(feature.location.start+1)}..{str(feature.location.end)})"
            else:
                #print(feature.join)
                pos = "complement(join(" + ','.join( f"{i.start+1}..{i.end}" for i in features[1].join.parts) + "))"
        
        if feature.type == "tRNA":
            gb_text += f'     tRNA            {pos}\n'
            gb_text += f'                     /product="{feature.name}"\n'
        elif feature.type == "rRNA":
            gb_text += f'     rRNA            {pos}\n'
            gb_text += f'                     /product="{product.get(feature.name, feature.name)}"\n'
        elif feature.type == "CDS":
            gb_text += f'     gene            {pos}\n'
            gb_text += f'                     /gene="{feature.name}"\n'
            gb_text += f'     CDS             {pos}\n'
            gb_text += f'                     /gene="{feature.name}"\n'
            gb_text += f'                     /codon_start=1\n'
            gb_text += f'                     /transl_table={table}\n'
            gb_text += f'                     /product="{product.get(feature.name, feature.name)}"\n'
            gb_text += f'                     {get_translation_string(feature, mtgenome=features[0].mtgenome, table=table)}\n'
        elif feature.type == "D-loop":
            gb_text += f'     D-loop          {pos}\n'
            gb_text += f'                     /note="Control Region"\n'

    if features[0].mtgenome == None:
        gb_text += f"CONTIG      join({organism}:1..16913)\n//"
    else:
        gb_text += f'ORIGIN\n{format_mtgenome(str(features[0].mtgenome).lower(), line_length=60, group_length=10)}\n//'

    if output == None:
        print(gb_text)
    else:
        with open(output, 'w') as f:
            f.write(gb_text)

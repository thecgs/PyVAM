#!/usr/bin/env python
# coding: utf-8

import re
import os
import time
import logging
from copy import copy
from Bio import SeqIO, Entrez, Data
from .config import CommonNamesDict, MTColors
from Bio.SeqFeature import CompoundLocation, ExactPosition, SimpleLocation, SeqFeature, Reference
from Bio.SeqRecord import SeqRecord
from functools import lru_cache

# Aliases may be removed when they are covered by the flexible matcher, but
# their canonical display names remain recognized results.
CANONICAL_GENE_NAMES = frozenset(CommonNamesDict.values())

class Feature:
    def __init__(self, name, location, type, color, join=None, mtgenome=None, accession=None, file=None, topology=None, partition=None, codon_start=1, locus_tags=(), original_type=None):
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
        self.codon_start = codon_start
        self.locus_tags = tuple(locus_tags)
        self.original_type = original_type
        
    def __repr__(self):
            return str(self.name)
    def __str__(self):
            return self.__repr__()


def alias_regex(alias):
    parts = re.split(r"[\s_/#-]+", alias.strip())
    return r"[\s_/#-]*".join(re.escape(part) for part in parts)

@lru_cache(maxsize=None)
def alias_pattern(alias):
    pattern = (
        r"(?:^|[\s_/#:,\-])"
        + alias_regex(alias)
        + r"(?=$|[\s_/#:,\-()])"
    )
    return re.compile(pattern)
    
def search_name(gene_name):
    gene_name = gene_name.upper()
    if gene_name in CommonNamesDict:
        return CommonNamesDict[gene_name]
    for alias, canonical in CommonNamesDict.items():
        if alias_pattern(alias).search(gene_name):
            return canonical
    return gene_name

def feature_key(feature):
    """Identify a biological annotation, including all its drawing parts."""
    location = feature.join if feature.join is not None else feature.location
    return (feature.name, feature.type, feature.locus_tags, feature.original_type,
            feature.codon_start, str(location))


def is_repeat(features, location, name=None, annotation=None):
    """Deduplicate exact annotations, not merely overlapping genomic spans."""
    feature_type = get_type(name) if name is not None else None
    for feature in features:
        if feature.type == "source" or (feature_type is not None and feature.type != feature_type):
            continue
        if annotation is not None:
            if feature.locus_tags != tuple(sorted(annotation.qualifiers.get("locus_tag", []))):
                continue
            if feature.original_type != annotation.type:
                continue
            original_location = feature.join if feature.join is not None else feature.location
            if original_location != annotation.location:
                continue
            if annotation.type == "CDS" and feature.codon_start != int(
                annotation.qualifiers.get("codon_start", ["1"])[0]
            ):
                continue
        if feature.location == location and (name is None or feature.name == name):
            return True
    return False


def _recognized_rna_name(annotation):
    """Prefer a known RNA gene name, then a product of the same RNA type."""
    for key in ("gene", "product"):
        for value in annotation.qualifiers.get(key, []):
            name = search_name(value)
            name = CommonNamesDict.get(name.upper(), name)
            if name in CANONICAL_GENE_NAMES and get_type(name) == annotation.type:
                return name
    return None


def _gene_has_annotation(gene, annotations):
    """Match gene placeholders to their CDS/RNA annotation by identity and span."""
    def names(annotation):
        values = annotation.qualifiers.get("gene") or annotation.qualifiers.get("product", [])
        return {CommonNamesDict.get(value.upper(), value)
                for value in (search_name(value) for value in values)}

    def contains(outer, inner):
        return all(any(a.start <= b.start and b.end <= a.end
                       and a.strand == b.strand and a.ref == b.ref
                       and a.ref_db == b.ref_db for a in outer.parts)
                   for b in inner.parts)

    if gene.location is None:
        return False
    gene_names = names(gene)
    gene_tags = set(gene.qualifiers.get("locus_tag", []))
    for annotation in annotations:
        if annotation.location is None or annotation.location.strand != gene.location.strand:
            continue
        annotation_tags = set(annotation.qualifiers.get("locus_tag", []))
        same_gene = (bool(gene_tags & annotation_tags) if gene_tags and annotation_tags
                     else bool(gene_names & names(annotation)))
        if same_gene and (contains(gene.location, annotation.location)
                          or contains(annotation.location, gene.location)):
            return True
    return False


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
    value = None
    if isinstance(start, str):
        for feature in features[1:]:
            if feature.name == start:
                value = feature.location.start
                break
                
    if value == None:
        logger.warning(
            "%s does not contain %s. Therefore, the rotation operation was not performed. Please check.",
            features[0].file,
            start,
        )
        return features
        
    if features[0].topology == "linear" and force_reoriented == False:
        logger.warning(
            "%s is linear topology. Reoriention is disabled by default, so no rotation operation was performed. To force reoriention, use the force_reoriented=True parameter.",
            features[0].file,
        )
        return features
        
    genome_length = len(features[0].mtgenome)

    def rotate_part(part):
        if part.end <= value:
            return [part + (genome_length - value)]
        if part.start >= value:
            return [part - value]
        # The new origin cuts this part. Keep the pieces in extraction order,
        # which is reversed for a feature on the negative strand.
        # Position addition preserves fuzzy boundary types; subtraction of
        # these int subclasses would discard the < or > annotation.
        parts = [
            SimpleLocation(part.start + (genome_length - value), genome_length,
                           strand=part.strand),
            SimpleLocation(0, part.end + (-value), strand=part.strand),
        ]
        return parts[::-1] if part.strand == -1 else parts

    source = copy(features[0])
    source.mtgenome = rotate_seq(source.mtgenome, value)
    rotated = []
    seen_joins = set()
    for feature in features[1:]:
        if feature.type == "Gap":
            marker = copy(feature)
            marker.location = feature.location - value
            rotated.append(marker)
            continue
        location = feature.join if feature.join is not None else feature.location
        if feature.join is not None:
            key = feature_key(feature)
            if key in seen_joins:
                continue
            seen_joins.add(key)
        parts = [piece for part in location.parts for piece in rotate_part(part)]
        joined = (CompoundLocation(parts, operator=getattr(location, "operator", "join"))
                  if len(parts) > 1 else None)
        for part in parts:
            shifted = copy(feature)
            shifted.location = part
            shifted.join = joined
            rotated.append(shifted)
    return [source] + sorted(rotated, key=lambda feature: feature.location.start)


def get_type(genename):
    if genename in ['ND1', 'ND2', 'ND3', 'ND4L', 'ND4', 'ND5', 'ND6', 'COX1', 'COX2', 'COX3', 'ATPase6', 'ATPase8', 'Cytb']:
        return "CDS"
    elif 'rRNA' in genename:
        return "rRNA"
    elif 'tRNA' in genename:
        return "tRNA"
    elif genename == "D-loop":
        return "D-loop"
    else:
        return "CDS"
        
def get_features(file, abbr=False, colors=None, isfilename2species=False, start=None,
                 force_reoriented=False, default_topology="circular"):
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
        default_topology: {str} topology to use when the GenBank record does
                                not declare one ("circular" by default).
    """
    if colors == None:
        colors = MTColors['MITOFISH']
    elif isinstance(colors, str):
        # Return the palette, rather than its name, when a theme is unknown.
        colors = MTColors.get(colors.upper(), MTColors['MITOFISH'])
    elif isinstance(colors, dict):
        pass
    else:
        raise TypeError("colors must be None, a theme name (str), or a color mapping (dict).")

    if not isinstance(default_topology, str) or default_topology.lower() not in {"circular", "linear"}:
        raise ValueError("default_topology must be either 'circular' or 'linear'.")
    default_topology = default_topology.lower()
            
    features = []
    
    # Own the input stream explicitly: SeqIO does not close caller-owned
    # handles, and parsing or annotation validation can exit before EOF.
    if os.path.exists(file):
        handle = open(file)
    else:
        handle = get_genbank_from_ncbi(file)

    with handle:
        for record in SeqIO.parse(handle, 'genbank'):
            if not record.features:
                raise ValueError(f"{file}: record {record.id} has no feature annotations.")
            if 'data_file_division' in record.annotations:
                partition = record.annotations['data_file_division']
            else:
                partition = "UNA"

            topology = record.annotations.get('topology')
            if not isinstance(topology, str) or topology.lower() not in {"circular", "linear"}:
                logger = logging.getLogger(__name__)
                logger.warning(
                    "%s does not declare a valid topology; assuming %s. Pass "
                    "default_topology to the calling API to override this.",
                    file, default_topology,
                )
                topology = default_topology
            topology = topology.lower()
            mtgenome = record.seq.upper()
            accession = record.id

            # Match the original annotations so CDS takes precedence regardless
            # of input order, gene name, or whether the CDS has multiple parts.
            genes = [feature for feature in record.features if feature.type == "gene"]
            annotations = []
            for annotation in record.features:
                if annotation.type == "D_loop":
                    annotation = copy(annotation)
                    annotation.type = "D-loop"
                if annotation.type == "CDS" and not annotation.qualifiers.get("gene"):
                    for gene in genes:
                        if _gene_has_annotation(gene, [annotation]):
                            # Keep CDS coordinates and qualifiers, but inherit its
                            # identity before discarding the matching gene entry.
                            annotation = copy(annotation)
                            annotation.qualifiers = dict(annotation.qualifiers)
                            for key in ("gene", "product", "note"):
                                if gene.qualifiers.get(key):
                                    if key == "gene" or not any(
                                        annotation.qualifiers.get(k) for k in ("product", "note")
                                    ):
                                        annotation.qualifiers[key] = list(gene.qualifiers[key])
                                    break
                            else:
                                if not any(annotation.qualifiers.get(k) for k in ("product", "note")):
                                    annotation.qualifiers["gene"] = list(gene.qualifiers["locus_tag"])
                            break
                annotations.append(annotation)
            typed_features = [feature for feature in annotations
                              if feature.type == "CDS" or (
                                  feature.type in ("tRNA", "rRNA")
                                  and _recognized_rna_name(feature) is not None)]
            for i in annotations:
                if i.type == "gene" and _gene_has_annotation(i, typed_features):
                    continue
                feature_start = len(features)
                if i.type == "CDS":
                    codon_start = int(i.qualifiers.get("codon_start", ["1"])[0])
                    if codon_start not in (1, 2, 3):
                        raise ValueError("CDS codon_start must be 1, 2, or 3.")
                # An explicit or inherited gene identifier takes precedence over
                # free-text descriptions, including for non-standard ORF names.
                if i.type != "source" and i.qualifiers.get("gene"):
                    gene_name = search_name(i.qualifiers["gene"][0])
                elif len(list(i.qualifiers.values())) != 0:
                    if "product" in i.qualifiers:
                        gene_name = i.qualifiers['product'][0]
                        gene_name =  search_name(gene_name)

                        #print("gene_product", gene_name)
                        if gene_name not in CANONICAL_GENE_NAMES and "gene" in i.qualifiers:
                            gene_name = i.qualifiers['gene'][0]
                            gene_name =  search_name(gene_name)

                            if gene_name not in CANONICAL_GENE_NAMES and "note" in i.qualifiers:
                                gene_name = i.qualifiers['note'][0]
                                gene_name =  search_name(gene_name)
                                if gene_name not in CANONICAL_GENE_NAMES and i.type.upper() in CommonNamesDict:
                                    gene_name =  search_name(i.type)

                    elif "gene" in i.qualifiers:
                        gene_name = i.qualifiers['gene'][0]
                        gene_name =  search_name(gene_name)
                        #print(i, gene_name)
                        if gene_name not in CANONICAL_GENE_NAMES and "note" in i.qualifiers:
                            #print(i, gene_name)
                            gene_name = i.qualifiers['note'][0]
                            gene_name =  search_name(gene_name)
                            if gene_name not in CANONICAL_GENE_NAMES and i.type.upper() in CommonNamesDict:
                                gene_name =  search_name(i.type)

                    elif "note" in i.qualifiers:
                        #print(i.qualifiers)
                        gene_name = i.qualifiers["note"][0]
                        gene_name =  search_name(gene_name)
                        #print(i, gene_name)
                        if gene_name not in CANONICAL_GENE_NAMES and i.type.upper() in CommonNamesDict:
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
                # RNA identifiers (e.g. MT-TF) need not be standard display names.
                # Recover a recognized product of the same RNA type without
                # changing the explicit identifiers of protein-coding ORFs.
                if i.type in ("tRNA", "rRNA"):
                    gene_name = _recognized_rna_name(i) or gene_name

                if i.type ==  "source":
                    if isfilename2species:
                        species_name = os.path.splitext(os.path.basename(file))[0]
                    else:
                        species_values = i.qualifiers.get("organism") or record.annotations.get("organism")
                        if not species_values:
                            raise ValueError(f"{file}: source feature is missing an organism name.")
                        species_name = species_values[0] if isinstance(species_values, (list, tuple)) else species_values

                    species_name = get_species_name(species_name, abbr=abbr)
                    features.append(Feature(name=species_name, location=i.location, type=i.type, color=colors.get('source', colors.get('Other genes', 'gray')),
                                            mtgenome=mtgenome, accession=accession, file=file, topology=topology, partition=partition))

                elif i.type in ['rRNA', 'tRNA', 'D_loop', 'D-loop']:
                    if gene_name in ['tRNA-His', 'tRNA-Pro', 'tRNA-Thr', 'tRNA-Trp', 'tRNA-Met', 'tRNA-Asp', 'tRNA-Ala', 'tRNA-Gln',
                                     'tRNA-Ile', 'tRNA-Arg', 'tRNA-Tyr', 'tRNA-Phe', 'tRNA-Lys', 'tRNA-Gly', 'tRNA-Asn', 'tRNA-Leu',
                                     'tRNA-Glu', 'tRNA-Val', 'tRNA-Cys', 'tRNA-Ser', '12S rRNA', '16S rRNA', "D-loop"]:
                        if isinstance(i.location, CompoundLocation):
                            for location in i.location.parts:
                                if not is_repeat(features=features, location=location, name=gene_name, annotation=i):
                                    features.append(Feature(name=gene_name, location=location, type=i.type, color=colors.get(gene_name, colors.get('Other genes', 'gray')),join=i.location))
                        else:
                            if not is_repeat(features=features, location=i.location, name=gene_name, annotation=i):
                                features.append(Feature(name=gene_name, location=i.location, type=i.type, color=colors.get(gene_name, colors.get('Other genes', 'gray'))))

                elif i.type in ['CDS', 'gene']:
                    if gene_name in ['ND1', 'ND2', 'ND3', 'ND4L', 'ND4', 'ND5', 'ND6', 'COX1', 'COX2', 'COX3', 'ATPase6', 'ATPase8', 'Cytb']:
                        if isinstance(i.location, CompoundLocation):
                            for location in i.location.parts:
                                if not is_repeat(features=features, location=location, name=gene_name, annotation=i):
                                    features.append(Feature(name=gene_name, location=location, type="CDS", color=colors.get(gene_name, colors.get('Other genes', 'gray')),join=i.location))
                        else:
                            if not is_repeat(features=features, location=i.location, name=gene_name, annotation=i):
                                features.append(Feature(name=gene_name, location=i.location, type="CDS", color=colors.get(gene_name, colors.get('Other genes', 'gray'))))

                    elif 'tRNA' in gene_name:
                        if isinstance(i.location, CompoundLocation):
                            for location in i.location.parts:
                                if not is_repeat(features=features, location=location, name=gene_name, annotation=i):
                                    features.append(Feature(name=gene_name, location=location, type="tRNA", color=colors.get(gene_name, colors.get('Other genes', 'gray')),join=i.location))
                        else:
                            if not is_repeat(features=features, location=i.location, name=gene_name, annotation=i):
                                features.append(Feature(name=gene_name, location=i.location, type="tRNA", color=colors.get(gene_name, colors.get('Other genes', 'gray'))))

                    elif gene_name in ['12S rRNA', '16S rRNA']:
                        if isinstance(i.location, CompoundLocation):
                            for location in i.location.parts:
                                if not is_repeat(features=features, location=location, name=gene_name, annotation=i):
                                    features.append(Feature(name=gene_name, location=location, type="rRNA", color=colors.get(gene_name, colors.get('Other genes', 'gray')),join=i.location))
                        else:
                            if not is_repeat(features=features, location=i.location, name=gene_name, annotation=i):
                                features.append(Feature(name=gene_name, location=i.location, type="rRNA", color=colors.get(gene_name, colors.get('Other genes', 'gray'))))
                    else: #ORF
                        if isinstance(i.location, CompoundLocation):
                            for location in i.location.parts:
                                if not is_repeat(features=features, location=location, name=gene_name, annotation=i):
                                    features.append(Feature(name=gene_name, location=location, type="CDS", color=colors.get(gene_name, colors.get('Other genes', 'gray')),join=i.location))
                        else:
                            if not is_repeat(features=features, location=i.location, name=gene_name, annotation=i):
                                features.append(Feature(name=gene_name, location=i.location, type="CDS", color=colors.get(gene_name, colors.get('Other genes', 'gray'))))


                elif i.type in ['misc_feature', 'repeat_region']:
                    if gene_name in ['tRNA-His', 'tRNA-Pro', 'tRNA-Thr', 'tRNA-Trp', 'tRNA-Met', 'tRNA-Asp', 'tRNA-Ala', 'tRNA-Gln',
                                     'tRNA-Ile', 'tRNA-Arg', 'tRNA-Tyr', 'tRNA-Phe', 'tRNA-Lys', 'tRNA-Gly', 'tRNA-Asn', 'tRNA-Leu',
                                     'tRNA-Glu', 'tRNA-Val', 'tRNA-Cys', 'tRNA-Ser']:
                        if isinstance(i.location, CompoundLocation):
                            for location in i.location.parts:
                                if not is_repeat(features=features, location=location, name=gene_name, annotation=i):
                                    features.append(Feature(name=gene_name, location=location, type="tRNA", color=colors.get(gene_name, colors.get('Other genes', 'gray')),join=i.location))
                        else:
                            if not is_repeat(features=features, location=i.location, name=gene_name, annotation=i):
                                features.append(Feature(name=gene_name, location=i.location, type="tRNA", color=colors.get(gene_name, colors.get('Other genes', 'gray'))))

                    elif gene_name in ['12S rRNA', '16S rRNA']: #, "D-loop"
                        if isinstance(i.location, CompoundLocation):
                            for location in i.location.parts:
                                if not is_repeat(features=features, location=location, name=gene_name, annotation=i):
                                    features.append(Feature(name=gene_name, location=location, type="rRNA", color=colors.get(gene_name, colors.get('Other genes', 'gray')),join=i.location))
                        else:
                            if not is_repeat(features=features, location=i.location, name=gene_name, annotation=i):
                                features.append(Feature(name=gene_name, location=i.location, type="rRNA", color=colors.get(gene_name, colors.get('Other genes', 'gray'))))

                    elif gene_name in ["D-loop"]:
                        if isinstance(i.location, CompoundLocation):
                            for location in i.location.parts:
                                if not is_repeat(features=features, location=location, name=gene_name, annotation=i):
                                    features.append(Feature(name=gene_name, location=location, type="D-loop", color=colors.get(gene_name, colors.get('Other genes', 'gray')),join=i.location))
                        else:
                            if not is_repeat(features=features, location=i.location, name=gene_name, annotation=i):
                                features.append(Feature(name=gene_name, location=i.location, type="D-loop", color=colors.get(gene_name, colors.get('Other genes', 'gray'))))


                    elif gene_name in ['ND1', 'ND2', 'ND3', 'ND4L', 'ND4', 'ND5', 'ND6', 'COX1', 'COX2', 'COX3', 'ATPase6', 'ATPase8', 'Cytb']:
                        if isinstance(i.location, CompoundLocation):
                            for location in i.location.parts:
                                if not is_repeat(features=features, location=location, name=gene_name, annotation=i):
                                    features.append(Feature(name=gene_name, location=location, type="CDS", color=colors.get(gene_name, colors.get('Other genes', 'gray')),join=i.location))
                        else:
                            if not is_repeat(features=features, location=i.location, name=gene_name, annotation=i):
                                features.append(Feature(name=gene_name, location=i.location, type="CDS", color=colors.get(gene_name, colors.get('Other genes', 'gray'))))
                else:
                    pass

                # Keep provenance and reading frame on each drawing part so
                # downstream grouping cannot merge independent annotations.
                for feature in features[feature_start:]:
                    feature.locus_tags = tuple(sorted(i.qualifiers.get("locus_tag", [])))
                    feature.original_type = i.type
                    feature.codon_start = codon_start if i.type == "CDS" else 1


    if not features:
        raise ValueError(f"{file}: no usable GenBank feature annotations were found.")
    source_index = next((index for index, feature in enumerate(features)
                         if feature.type == "source"), None)
    if source_index is None:
        raise ValueError(f"{file}: a usable source feature is required before gene annotations.")
    if source_index:
        features = [features[source_index]] + features[:source_index] + features[source_index + 1:]

    
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
    
def tidy_genbank(file, output=None, isfilename2species=False, start=None, table=2,
                 force_reoriented=False, partition="inherit", default_topology="circular"):
    """
    Descripton:
        Use PyVAM's powerful GenBank parser to reorganize the GenBank
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
        default_topology: {str} topology to use when the input record does not
                                declare one ("circular" or "linear").
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
    
    def get_translation_string(feature, mtgenome, table):
        # ``feature.location`` is the first part retained for drawing a
        # compound feature.  Use its original compound location when
        # translating so the exported protein covers every CDS segment.
        location = feature.join if feature.join is not None else feature.location
        CDS = location.extract(mtgenome)[feature.codon_start - 1:]
        CDS = CDS[:len(CDS) - (len(CDS) % 3)]
        pep = str(CDS.translate(table=table))
        # CompoundLocation parts are in extraction order, even across the
        # origin. The negative strand starts at the first part's end.
        first_part = location.parts[0]
        five_prime = first_part.end if first_part.strand == -1 else first_part.start
        complete_start = (first_part.strand in (1, -1)
                          and type(five_prime) is ExactPosition)
        if (feature.codon_start == 1 and complete_start
                and str(CDS[0:3]) in Data.CodonTable.unambiguous_dna_by_id[table].start_codons):
            pep = "M" + pep[1:]
        if pep.endswith("*"):
            pep = pep[:-1]
        return pep

    features = get_features(file, isfilename2species=isfilename2species, start=start,
                            force_reoriented=force_reoriented,
                            default_topology=default_topology)
    
    features_tmp = []
    tmp = []
    for feature in features:
        if feature.join == None:
            features_tmp.append(feature)
        else:
            key = feature_key(feature)
            if key not in tmp:
                features_tmp.append(feature)
                tmp.append(key)
    
    features = features_tmp
    organism = features[0].name
    genome_len = len(features[0].mtgenome)
    
    #print(features)
    
    if partition == "inherit":
        partition = features[0].partition
        
    record = SeqRecord(features[0].mtgenome, id=".",
                       name=organism.replace(" ", "_"), description=".")
    reference = Reference()
    reference.location = [SimpleLocation(0, genome_len)]
    reference.authors = "Chen, G."
    reference.title = "PyVAM: A Python package for visualizing animal mitochondrial."
    reference.journal = "Unpublished"
    record.annotations = {
        "molecule_type": "DNA",
        "topology": features[0].topology,
        "data_file_division": partition,
        "date": time.strftime("%d-%b-%Y", time.localtime()).upper(),
        "source": f"mitochondrion {organism}",
        "organism": organism,
        "taxonomy": ["Unclassified"],
        "references": [reference],
    }
    record.features = [SeqFeature(
        SimpleLocation(0, genome_len, strand=1), type="source",
        qualifiers={"organism": [organism], "organelle": ["mitochondrion"],
                    "mol_type": ["genomic DNA"]},
    )]

    for feature in features[1:]:
        location = feature.join if feature.join is not None else feature.location
        identity = {"locus_tag": list(feature.locus_tags)} if feature.locus_tags else {}
        if feature.type == "tRNA":
            qualifiers = {"product": [feature.name]}
        elif feature.type == "rRNA":
            qualifiers = {"product": [product.get(feature.name, feature.name)]}
        elif feature.type == "CDS":
            record.features.append(SeqFeature(
                location, type="gene", qualifiers={"gene": [feature.name], **identity},
            ))
            qualifiers = {
                "gene": [feature.name],
                "codon_start": [str(feature.codon_start)],
                "transl_table": [str(table)],
                "product": [product.get(feature.name, feature.name)],
                "translation": [get_translation_string(feature, features[0].mtgenome, table)],
            }
        elif feature.type == "D-loop":
            qualifiers = {"note": ["Control Region"]}
        else:
            continue
        qualifiers.update(identity)
        record.features.append(SeqFeature(location, type=feature.type, qualifiers=qualifiers))

    if output == None:
        print(record.format("genbank"), end="")
    else:
        SeqIO.write(record, output, "genbank")

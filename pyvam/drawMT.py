#!/usr/bin/env python
# coding: utf-8

import math, re
import logging
import matplotlib.pyplot as plt
from matplotlib import collections
from matplotlib.patches import Patch, FancyArrow
from matplotlib.lines import Line2D
from .parserGB import get_features
from .config import MTColors_legends, FullName2AbbrName

def rotation_text(theta):
    d =  (theta * 180 / math.pi)
    if d <180:
        return 90 - d 
    else:
        return  90 -d + 180

def get_GC(seq):
    try:
        return round((seq.count("G") + seq.count("C")) / (seq.count("A") + seq.count("T") + seq.count("G") + seq.count("C"))*100, 2)
    except ZeroDivisionError:
        return 0 

def get_GC_bar_param(features, bin=50, step=50):
    x = []
    y = []
    w = []
    scale_factor = (2 * math.pi) / (len(features[0].mtgenome) + 1)
    for i in range(0, len(features[0].mtgenome), step):
        if (i + bin) < len(features[0].mtgenome):
            start = i + 1
            end = i + bin + 1
            theta = ((end - start) / 2 + start) * scale_factor
            width = (end - start) * scale_factor
            GC = get_GC(features[0].mtgenome[i: i+bin])
            y.append(GC/100)
            x.append(theta)
            w.append(width)
        else:
            start = i + 1
            end = len(features[0].mtgenome) + 1
            theta = ((end - start) / 2 + start) * scale_factor
            width = (end - start) * scale_factor
            GC = get_GC(features[0].mtgenome[i: len(features[0].mtgenome)])
            y.append(GC/100)
            x.append(theta)
            w.append(width)
    return x, y, w
    
def stat_features(features):
    tRNA_count = 0
    rRNA_count = 0
    PCG_count = 0
    for feature in features:
        if re.search('tRNA', feature.name) and feature.join==None:
            tRNA_count += 1
        elif re.search('rRNA', feature.name) and feature.join==None:
            rRNA_count += 1
        elif feature.type == 'CDS' and feature.join==None:
        #elif re.search('COX|ATP|Cytb|ND', feature.name) and feature.join==None:
            PCG_count += 1

    tmp = []
    features_join = []
    for feature in features:
        if feature.join not in tmp:
            tmp.append(feature.join)
            features_join.append(feature)

    for feature in features_join:
        if re.search('tRNA', feature.name):
            tRNA_count += 1
        elif re.search('rRNA', feature.name):
            rRNA_count += 1
        elif feature.type == 'CDS':
        #elif re.search('COX|ATP|Cytb|ND', feature.name):
            PCG_count += 1
    return tRNA_count, rRNA_count, PCG_count

def draw_circos_MT(file,
                   output=None,
                   abbr=False,
                   isfilename2species=False,
                   colors="mitofish",
                   radius=25,
                   show_gene_label = True,
                   gene_label_size=7,
                   gene_label_inner=False,
                   show_info=True,
                   info_fontsize=15,
                   show_legend=True,
                   legend_size=6,
                   legend_postion=(1, -0.1),
                   show_GC_circos=True,
                   GC_circos_height = 0.3,
                   GC_circos_color = 'grey',
                   GC_circos_bin = 50,
                   GC_circos_step = 50,
                   start=None,
                   axes=None,
                   direction=-1,
                   figsize=(10,10),
                   tidyname=False,
                   add_id = False,
                   dpi = 300,
                   remove_NCR = False,
                  ):
    """
    Descripton:
        Draw a mitochondrial circos.
        
    Parameters：
        file: {str} a genbankfile or NCBI accession ID.
        output: {str} a path of figure save.
        abbr: {bool} whether to abbreviate species names.
        isfilename2species: {bool} whether filename convert to species.
        colors: {str, dict} themes such as, Chen, Tan, ogdraw, mitofish,
                            mitofish1, mitoz,  gggenes, chloroplot, grey, igv.
        radius: {int} radius of the circle.
        show_gene_label: {bool} show gene label.
        gene_label_size: {int} gene label fontsize.
        gene_label_inner: {bool} whether the control gene label is inner.
        show_info: {bool} show species name, GC context, genome length, and gene count information.
        info_fontsize: {int} information fontsize.
        show_legend: {bool} show fig legend.
        legend_size: {int} legend size.
        legend_postion: {tuple} such as (x, y), x=0-1, y=0-1.
        show_GC_circos: {bool} show GC circos.
        GC_circos_height {0-1} height of GC_circos.
        GC_circos_color {color} color of GC_circos.
        GC_circos_step {int} step of GC_circos (bp).
        GC_circos_bin {int} bin of GC_circos (bp).
        start: {None, str} initial feature, such as, ND1, ND2, ND3, ND4, ND4L, ND5, ND6,
                     COX1, COX2, COX3, ATPase6, ATPase8, Cytb, tRNA-His, tRNA-Pro,
                     tRNA-Thr, tRNA-Trp, tRNA-Met, tRNA-Asp, tRNA-Ala, tRNA-Gln,
                     tRNA-Ile, tRNA-Arg, tRNA-Tyr, tRNA-Phe, tRNA-Lys, tRNA-Gly,
                     tRNA-Asn, tRNA-Leu, tRNA-Glu, tRNA-Val, tRNA-Cys, tRNA-Ser,
                     12S rRNA, 16S rRNA, D-loop. default=None.
        axes: {matplotlib.projections.polar.PolarAxes}
        figsize: {tuple} fig of size.
        direction: {1, -1} clockwise: -1, anticlockwise: 1
        tidyname: {bool} tidy gene name.
        add_id: {bool} Species add to accession id from NCBI.
        dpi: {int} dpi value. the resolution in dots per inch.
        remove_NCR {bool}: Do not display D-loop; this may be useful for comparing a mix of genbank with and without NCR annotation.
    """
    add_other_genes = False
    if axes==None:
        fig, ax = plt.subplots(1,1, subplot_kw={'projection':'polar'}, figsize=figsize)
    else:
        ax = axes
    
    ax.set_theta_zero_location('N')
    #ax.set_theta_offset(offset)
    ax.set_theta_direction(direction)
    ax.grid(False)
    ax.set_frame_on(False)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    
    features = get_features(file, abbr=abbr, isfilename2species=isfilename2species, colors=colors, start=start)
    if remove_NCR:
        features = [feature for feature in features if feature.name != "D-loop"]
    add_other_genes = is_othergenes(features)
    try:
        #str(features[0].mtgenome)  # UndefinedSequenceError
        GC = get_GC(features[0].mtgenome)
        genome_length = len(features[0].mtgenome)
        tRNA_count, rRNA_count, PCG_count = stat_features(features)
        if add_id:
            info = features[0].accession +"; " + features[0].name+f"\n{genome_length} bp\nGC: {GC}%\n{PCG_count} PCGs; {rRNA_count} rRNAs; {tRNA_count} tRNAs"
        else:
            info = features[0].name+f"\n{genome_length} bp\nGC: {GC}%\n{PCG_count} PCGs; {rRNA_count} rRNAs; {tRNA_count} tRNAs"
    except:
        #show_info = False
        show_GC_circos = False
        logger = logging.getLogger(__name__) 
        logger.setLevel(logging.DEBUG)
        logger.warning(features[0].file+" Sequence content is undefined.")
        tRNA_count, rRNA_count, PCG_count = stat_features(features)
        if add_id:
            info = features[0].accession +"; " + features[0].name+f"\n{PCG_count} PCGs; {rRNA_count} rRNAs; {tRNA_count} tRNAs"
        else:
            info = features[0].name+f"\n{PCG_count} PCGs; {rRNA_count} rRNAs; {tRNA_count} tRNAs"
        
    if show_info:            
        #GC = get_GC(features[0].mtgenome)
        #genome_length = len(features[0].mtgenome)
        ax.text(0.5, 0.5, s=info,
                fontsize=info_fontsize, style='italic', va='center', ha='center')
                
    if tidyname:
        for i, f in enumerate(features):
            features[i].name = FullName2AbbrName.get(f.name, f.name)
    
    
    if features[0].topology == 'linear':
        logger = logging.getLogger(__name__) 
        logger.setLevel(logging.DEBUG)
        logger.warning(features[0].file+" is linear topology. Therefore, the rotation operation was not performed.")
        scale_factor = (1.9 * math.pi) / (int(features[0].location.end) + 1)
        ax.set_thetamin(0)
        ax.set_thetamax(342)        
    else:
        scale_factor = (2 * math.pi) / (int(features[0].location.end) + 1)
            
    for feature in features[1:]:
        start = int(feature.location.start) + 1
        end = int(feature.location.end) + 1
        theta = ((end - start) / 2 + start) * scale_factor
        width = (end - start) * scale_factor
        
        if gene_label_inner==False:
            if feature.location.strand == 1:
                ax.bar(theta, 1, width=width, bottom=radius-5, color=feature.color, label=feature.name, linewidth=0.5, edgecolor='black')
                if show_gene_label:
                    ax.annotate(feature.name,  xy=(theta, radius-4), xytext=(theta, radius-1), arrowprops=dict(arrowstyle="-", connectionstyle="arc3"), 
                                horizontalalignment='center', verticalalignment='center', fontsize=gene_label_size, rotation=rotation_text(theta))
            elif feature.location.strand == -1:
                ax.bar(theta, 1, width=width, bottom=radius-6, color=feature.color, label=feature.name, linewidth=0.5, edgecolor='black')
                if show_gene_label:
                    ax.annotate(feature.name,  xy=(theta, radius-5), xytext=(theta, radius), arrowprops=dict(arrowstyle="-", connectionstyle="arc3"), 
                                horizontalalignment='center', verticalalignment='center', fontsize=gene_label_size, rotation=rotation_text(theta))
            ax.plot([i* math.pi/180 for i in range(0, 361)], [radius-5]*361, c='black')
        
        else:
            if feature.location.strand == 1:
                ax.bar(theta, 1, width=width, bottom=radius, color=feature.color, label=feature.name, linewidth=0.5, edgecolor='black')
                if show_gene_label:
                    ax.annotate(feature.name,  xy=(theta, radius), xytext=(theta, radius-4), arrowprops=dict(arrowstyle="-", connectionstyle="arc3"), 
                                horizontalalignment='center', verticalalignment='center', fontsize=gene_label_size, rotation=rotation_text(theta))
            else:
                ax.bar(theta, 1, width=width, bottom=radius-1, color=feature.color, label=feature.name, linewidth=0.5, edgecolor='black')
                if show_gene_label:
                    ax.annotate(feature.name,  xy=(theta, radius-1), xytext=(theta, radius-5), arrowprops=dict(arrowstyle="-", connectionstyle="arc3"), 
                                horizontalalignment='center', verticalalignment='center', fontsize=gene_label_size, rotation=rotation_text(theta))
                
            ax.plot([i* math.pi/180 for i in range(0, 361)], [radius]*361, c='black')
            
    if colors == None:
        legend_elements = MTColors_legends['MITOFISH']
        ncol = 1
    elif isinstance(colors, str):
        legend_elements = MTColors_legends.get(colors.upper(), 'MITOFISH')
        if colors.upper() == "IGV":
            ncol = 4
        else:
            ncol = 1
    elif isinstance(colors, dict):
        legend_elements = [Patch(facecolor=colors[i], edgecolor='black', label=i) for i in colors if i!="source"]
        ncol = 4
    
    if add_other_genes == False:
        legend_elements = legend_elements[:-1]
    
    if show_legend:
        ax.legend(handles=legend_elements,
                  loc='center', 
                  bbox_to_anchor=legend_postion, 
                  ncol=ncol,
                  #frameon=False,
                  #shadow=False, 
                  prop={'size': legend_size},
                  title='')
    
    if show_GC_circos:
        x,y,w = get_GC_bar_param(features, step=GC_circos_step, bin=GC_circos_bin)
        ax.bar(x, [i*10 *GC_circos_height for i in y], width=w, bottom=radius-10, color=GC_circos_color, linewidth=0, edgecolor='black')
    
    if output!=None:
        fig.savefig(output, bbox_inches='tight', dpi=dpi)
    
    if axes==None:
        return fig, ax
 
def _get_mt_rect_ax(ax, features, show_info=True, info_fontsize=10, show_gene_label=True, gene_label_size=5, add_id=False, force_reoriented=True, species_label_size=10, species_label_color='black'):
    height = 0.2     # gene width
    linewidth = 0.03 # genome width
    _staus_brake_tmp = False
    _staus_brake = False
    box_width = 200
    
    for index, feature in enumerate(features):
        start = int(feature.location.start)
        end = int(feature.location.end)
        if _staus_brake_tmp:
            start += box_width
            end += box_width
            
        if feature.type == 'source':
            rect = plt.Rectangle(xy=(0, 0.45), width=end, height=linewidth,  color=feature.color)
            species = feature.name
            genome_length = int(feature.location.end)
            
        elif feature.location.strand == 1:
            rect = plt.Rectangle(xy=(start, 0.45+linewidth), 
                                 width=end-start, 
                                 height=height, 
                                 facecolor=feature.color,
                                 edgecolor='black')
            if show_gene_label:
                x = (end-start)/2 + start
                ax.text(x, y=0.45+height+linewidth+0.05, s=feature.name, rotation = 'vertical', size=gene_label_size, ha='center')
            
        elif feature.location.strand == -1:
            rect = plt.Rectangle(xy=(start, 0.25), 
                                 width=end-start,
                                 height=height, 
                                 facecolor=feature.color,
                                 edgecolor='black')
            if show_gene_label:
                x = (end-start)/2 + start
                ax.text(x, y=0.45-height-linewidth+0.005, s=feature.name, rotation = 'vertical', size=gene_label_size, va='top', ha='center')
        else:
            _staus_brake = True
            _staus_brake_tmp = True
            x = end
            line = collections.LineCollection([[(x+(box_width/6), 0.25), (x+(box_width/2), 0.65)],
                                               [(x+(box_width/2), 0.25), (x+(box_width*5/6), 0.65)]],
                                              color=["black", "black"])
            ax.add_collection(line)
            if features[-1].name != "Gap":
                rect1 = plt.Rectangle(xy=(genome_length, 0.45), width=box_width, height=linewidth, color='black')
                rect2 = plt.Rectangle(xy=(x+7, 0.42), width=box_width-16, height=linewidth+0.05, edgecolor=None, facecolor='white', linewidth=0)
                ax.add_patch(rect1)
                ax.add_patch(rect2)
            
        ax.add_patch(rect)
        
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    if add_id:
        ax.text(-100, 0.5, s=species+" "+features[0].accession, ha='right', style='italic', fontsize=species_label_size, color=species_label_color)
    else:
        ax.text(-100, 0.5, s=species, ha='right', style='italic', fontsize=species_label_size, color=species_label_color)
        
    if show_info:
        if _staus_brake_tmp:
            ax.text(genome_length+100+box_width, 0.55, s="+", ha='center', weight='bold', fontsize=info_fontsize)
            ax.text(genome_length+100+box_width, 0.25, s="-", ha='center', weight='bold', fontsize=info_fontsize)
            ax.text(genome_length+200+box_width, 0.47, s=str(genome_length)+'bp', ha='left', fontsize=info_fontsize)
        else:
            ax.text(genome_length+100, 0.55, s="+", ha='center', weight='bold', fontsize=info_fontsize)
            ax.text(genome_length+100, 0.25, s="-", ha='center', weight='bold', fontsize=info_fontsize)
            ax.text(genome_length+200, 0.47, s=str(genome_length)+'bp', ha='left', fontsize=info_fontsize)
    return _staus_brake

def is_othergenes(features):
    res = False
    for feature in features[1:]:
        if feature.name not in ['ND1', 'ND2', 'ND3', 'ND4L', 'ND4', 'ND5', 'ND6',
     'COX1', 'COX2', 'COX3', 'ATPase6', 'ATPase8', 'Cytb', 'tRNA-His',
     'tRNA-Pro', 'tRNA-Thr', 'tRNA-Trp', 'tRNA-Met', 'tRNA-Asp',
     'tRNA-Ala', 'tRNA-Gln', 'tRNA-Ile', 'tRNA-Arg', 'tRNA-Tyr', 
     'tRNA-Phe', 'tRNA-Lys', 'tRNA-Gly', 'tRNA-Asn', 'tRNA-Leu',
     'tRNA-Glu', 'tRNA-Val', 'tRNA-Cys', 'tRNA-Ser', '12S rRNA',
     '16S rRNA', 'D-loop', 'Gap']:
            #print(feature.name)
            res = True
    return res

def draw_linear_MT(files,
                   output=None,
                   abbr=False,
                   isfilename2species=False,
                   colors="mitofish",
                   show_gene_label = True,
                   gene_label_size=5,
                   gene_label_color='black',
                   species_label_size=10,
                   species_label_color='black',
                   show_info=True,
                   info_fontsize=10,
                   show_legend=True,
                   legend_size=10,
                   legend_postion=(1, 0.5),
                   start=None,
                   show_xaxis=True,
                   xaxisfontsize=10, 
                   hspace=0.5,
                   subplot_height_cex=1.5,
                   tidyname=False,
                   add_id = False,
                   dpi = 300,
                   force_reoriented=True,
                   remove_NCR = False,
                  ):
    
    """
    Descripton:
        The order of mitochondrial genes, arranged in proportion to their size.
    
    Parameters：
        files: {str, list, tuple} one or more genbankfile or NCBI accession ID.
        output: {str} a path of fig save.
        abbr: {bool} whether to abbreviate species names.
        isfilename2species: {bool} whether filename convert to species.
        colors: {str, dict} themes such as, Chen, Tan, ogdraw, mitofish,
                            mitofish1, mitoz,  gggenes, chloroplot, grey, igv.

        show_gene_label: {bool} show gene label.
        gene_label_size: {int} gene label size.
        gene_label_color: {str} gene label color.
        species_label_size: {int} species name size.
        species_label_color: {str} species name color.
        show_info: {bool} show genome length, and strand.
        info_fontsize: {int} information fontsize.
        show_legend: {bool} show fig legend.
        legend_size: {int} legend size.
        legend_postion: {tuple} such as (x, y), x=0-1, y=0-1.
        start: {None, str} initial feature, such as, ND1, ND2, ND3, ND4, ND4L, ND5, ND6,
                     COX1, COX2, COX3, ATPase6, ATPase8, Cytb, tRNA-His, tRNA-Pro,
                     tRNA-Thr, tRNA-Trp, tRNA-Met, tRNA-Asp, tRNA-Ala, tRNA-Gln,
                     tRNA-Ile, tRNA-Arg, tRNA-Tyr, tRNA-Phe, tRNA-Lys, tRNA-Gly,
                     tRNA-Asn, tRNA-Leu, tRNA-Glu, tRNA-Val, tRNA-Cys, tRNA-Ser,
                     12S rRNA, 16S rRNA, D-loop. default=None.
        show_xaxis: {bool} show x-axis.
        xaxisfontsize: {int} x-axis trick label fontsize.
        hspace: {0-1} spacing between subplots.
        subplot_height_cex: {float} a factors that control subplot scaling.
        tidyname: {bool} tidy gene name.
        add_id: {bool} Species add to accession id from NCBI.
        dpi: {int} dpi value. the resolution in dots per inch.
        force_reoriented: {bool} force-reoriendted linear mtgenome.
        remove_NCR {bool}: Do not display D-loop; this may be useful for comparing a mix of genbank with and without NCR annotation.
    """    
    
    add_other_genes = False
    
    if not (isinstance(files, list)) and (not isinstance(files, tuple)):
        files = [files]
    
    genome_max_length = 0
    genomes = []
    for file in files:
        features = get_features(file, abbr=abbr, isfilename2species=isfilename2species, colors=colors, start=start, force_reoriented=force_reoriented)
        if remove_NCR:
            features = [feature for feature in features if feature.name != "D-loop"]
        
        add_other_genes = is_othergenes(features)
        if tidyname:
            for i, f in enumerate(features):
                features[i].name = FullName2AbbrName.get(f.name, f.name)
        genomes.append(features)
        if genome_max_length < len(features[0].location):
            genome_max_length = len(features[0].location)
            
    if show_xaxis:
        fig, ax= plt.subplots(len(genomes)+1, 1, figsize=(20, len(genomes)*subplot_height_cex), sharex=True, sharey=True)
    else:
        fig, ax= plt.subplots(len(genomes), 1, figsize=(20, len(genomes)*subplot_height_cex), sharex=True, sharey=True)
        
    plt.subplots_adjust(hspace=hspace)
    
    if show_xaxis == False and len(genomes) == 1:
        for i, features in enumerate(genomes):
            _staus_brake = _get_mt_rect_ax(ax=ax, features=features, show_info=show_info,
                           info_fontsize=info_fontsize, show_gene_label=show_gene_label, 
                           gene_label_size=gene_label_size, add_id=add_id,force_reoriented=force_reoriented, species_label_size=species_label_size, species_label_color=species_label_color)
    else:
        for i, features in enumerate(genomes):
            _staus_brake = _get_mt_rect_ax(ax=ax[i], features=features, show_info=show_info,
                            info_fontsize=info_fontsize, show_gene_label=show_gene_label,
                            gene_label_size=gene_label_size, add_id=add_id,force_reoriented=force_reoriented, species_label_size=species_label_size, species_label_color=species_label_color)
        
    if show_xaxis:
        ax[-1].get_yaxis().set_visible(False)
        ax[-1].spines['left'].set_visible(False)
        ax[-1].spines['right'].set_visible(False)
        ax[-1].spines['bottom'].set_visible(False)
        ax[-1].set_xticks([i for i in range(0 , genome_max_length+2000, 1000)], labels=[i for i in range(0 , genome_max_length+2000, 1000)], **{"fontsize":xaxisfontsize})
        ax[-1].tick_params(top=True, labeltop=True, bottom=False, labelbottom=False, direction='in', pad=-20)
    
    if colors == None:
        legend_elements = MTColors_legends['MITOFISH']
    elif isinstance(colors, str):
        legend_elements = MTColors_legends.get(colors.upper(), 'MITOFISH')

    if show_legend:
        if _staus_brake:
            legend_elements = [Line2D([0], [0], color='black',marker='>', markersize=legend_size, 
                                      linestyle='-', linewidth=0, markerfacecolor='white', 
                                      markeredgecolor='white', label="// Break"
                                     )
                              ] + legend_elements
            
        if add_other_genes == False:
            legend_elements = legend_elements[:-1]
            
        ax[-1].legend(handles=legend_elements,
                      loc='upper right', 
                      bbox_to_anchor=legend_postion, 
                      ncol=12,
                      #frameon=False,
                      #shadow=False,
                      prop={'size': legend_size},
                      title='')
    
    if output!=None:
        fig.savefig(output, bbox_inches='tight', dpi=dpi)    
    return fig, ax
    
def get_box_param(genename):
    head_length=0.0025
    if genename in ["F","V","L","I","Q","M","W","A","N","C","Y","S","D","K","G","R","H","S","L","E","T","P"]:
        tail_length = 0.005
        text_offset = tail_length/2.5
        
    elif genename in ["ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"]:
        tail_length = 0.008
        if genename=="ND4L":
            text_offset = tail_length/8
        else:
            text_offset = tail_length/6
        
    elif genename in ["12S", "16S"]:
        tail_length = 0.005
        text_offset = tail_length/6
    elif genename in ['COX1','COX2','COX3']:
        tail_length = 0.008
        text_offset = tail_length/8
    elif genename in ["ATP6", "ATP8"]:
        tail_length = 0.008
        text_offset = tail_length/6
    elif genename == "NCR":
        tail_length = 0.007
        text_offset = tail_length/6
    elif genename == "Cytb":
        tail_length = 0.008
        text_offset = tail_length/5
    elif genename == 'Gap':
        tail_length = 0.005
        text_offset = tail_length/6
    else:
        tail_length = 0.015
        text_offset = tail_length/4
    box_width = head_length + tail_length
    return head_length, tail_length, box_width, text_offset
    
    
def remove_join(features):
    tmp = []
    res = []
    for feature in features:
        if feature.join!=None:
            if feature.join not in tmp:
                res.append(feature)
                tmp.append(feature.join)
        else:
            res.append(feature)
    return res
    
def draw_linear_MT_nonproportional(files,
                                   output=None,
                                   abbr=False,
                                   isfilename2species=False,
                                   colors="mitofish",
                                   gene_label_size=9, 
                                   gene_label_color="black",
                                   show_legend=True,
                                   legend_size=8,
                                   legend_postion=(1, 0),
                                   start=None,
                                   height=0.6,
                                   species_label_size=10,
                                   species_label_color="black",
                                   axes=None,
                                   add_id = False,
                                   dpi = 300,
                                   force_reoriented=False,
                                   remove_NCR = False,
                                  ):
    """
    Descripton:
        The order of mitochondrial genes is not depicted in proportion to their gene size.
    
    Parameters：
        files: {str, list, tuple} one or more genbankfile or NCBI accession ID.
        output: {str} a path of fig save.
        abbr: {bool} whether to abbreviate species names.
        isfilename2species: {bool} whether filename convert to species.
        colors: {str, dict} themes such as, Chen, Tan, ogdraw, mitofish,
                            mitofish1, mitoz,  gggenes, chloroplot, grey, igv.
        gene_label_size: {int} gene label size.
        gene_label_color: {str} gene label color.
        species_label_size: {int} species name size.
        species_label_color: {str} species name color.
        show_legend: {bool} show fig legend.
        legend_size: {int} legend size.
        legend_postion: {tuple} such as (x, y), x=0-1, y=0-1.
        start: {None, str} initial feature, such as, ND1, ND2, ND3, ND4, ND4L, ND5, ND6,
                     COX1, COX2, COX3, ATPase6, ATPase8, Cytb, tRNA-His, tRNA-Pro,
                     tRNA-Thr, tRNA-Trp, tRNA-Met, tRNA-Asp, tRNA-Ala, tRNA-Gln,
                     tRNA-Ile, tRNA-Arg, tRNA-Tyr, tRNA-Phe, tRNA-Lys, tRNA-Gly,
                     tRNA-Asn, tRNA-Leu, tRNA-Glu, tRNA-Val, tRNA-Cys, tRNA-Ser,
                     12S rRNA, 16S rRNA, D-loop. default=None.
        height: {0-1} box height.
        axes: {matplotlib.projections.polar.PolarAxes}
        add_id: {bool} Species add to accession id from NCBI.
        dpi: {int} dpi value. the resolution in dots per inch.
        force_reoriented: {bool} force-reoriendted linear mtgenome.
        remove_NCR {bool}: Do not display D-loop; this may be useful for comparing a mix of genbank with and without NCR annotation.
    """
    add_other_genes = False
    if not (isinstance(files, list)) and (not isinstance(files, tuple)):
        files = [files]
    
    genomes = []
    for file in reversed(files):
        features = get_features(file, abbr=abbr, isfilename2species=isfilename2species, colors=colors, start=start, force_reoriented=force_reoriented)
        if remove_NCR:
            features = [feature for feature in features if feature.name != "D-loop"]
        features = remove_join(features)
        genomes.append(features)
        add_other_genes = is_othergenes(features)
    
    if axes == None:
        fig, ax= plt.subplots(1, 1, figsize=(20, len(genomes)/2))
    else:
        ax = axes
    ax.set_ylim(0,len(genomes))
    
    xlimmax = 0
    ys = []
    ss = []
    _staus_brake = False
    for y, features in enumerate(genomes):
        for i, f in enumerate(features):
            features[i].name = FullName2AbbrName.get(f.name, f.name)
            
        if features[1].location.strand == -1:
            x = get_box_param(features[1].name)[2]+0
        else:
            x=0
            
        ys.append(y+0.5)
        if add_id:
            ss.append(features[0].name+ " " + features[0].accession)
        else:
            ss.append(features[0].name)
            
        for i, feature in enumerate(features[1:]):
            head_length, tail_length, box_width, text_offset= get_box_param(feature.name)
            
            if feature.location.strand == 1:
                arrow = FancyArrow(x, y+0.5, tail_length, 0, width=height, fc=feature.color, ec='black', head_width=height, head_length=head_length)
                ax.text(x+text_offset, y+0.5, feature.name, size=gene_label_size, color=gene_label_color, ha="left", va="center",clip_on=False, clip_path=None)
                
                if i+1 < len(features[1:]):
                    if features[1:][i+1].location.strand == -1:
                        x += get_box_param(features[1:][i+1].name)[2] + get_box_param(features[1:][i].name)[2]
                    else:
                        x += get_box_param(features[1:][i].name)[2]
            elif feature.location.strand == -1:
                arrow = FancyArrow(x, y+0.5, -tail_length, 0, width=height, fc=feature.color, ec='black', head_width=height, head_length=head_length)
                ax.text(x-text_offset, y+0.5, feature.name, size=gene_label_size, color=gene_label_color, ha="right", va="center", clip_on=False, clip_path=None)
                if i+1 < len(features[1:]):
                    if features[1:][i+1].location.strand == -1:
                        x +=  get_box_param(features[1:][i+1].name)[2]
            else:
                _staus_brake = True
                line = collections.LineCollection([[(x+box_width/6, y+0.2), (x+box_width/2, y+0.8)],
                                                   [(x+box_width/2, y+0.2), (x+box_width*5/6, y+0.8)]],
                                                  color=["black", "black"])
                ax.add_collection(line)
                
                if i+1 < len(features[1:]):
                    if features[1:][i+1].location.strand == -1:
                        x += get_box_param(features[1:][i+1].name)[2] + get_box_param(features[1:][i].name)[2]
                    else:
                        x += get_box_param(features[1:][i+1].name)[2]
            
            ax.add_patch(arrow)
            if xlimmax < x:
                xlimmax = x
                
    ax.set_xlim(0, xlimmax + xlimmax*0.05)           
    ax.get_xaxis().set_visible(False)
    #ax.get_yaxis().set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.yaxis.set_ticks_position('none')
    ax.set_yticks(ys)
    ax.set_yticklabels(ss, fontdict={'fontsize':species_label_size, 'fontstyle':'italic', 'color':species_label_color})
    
    if show_legend:
        if _staus_brake:
            legend_elements = [Line2D([0], [0], color='black',marker='>', markersize=legend_size, 
                                      linestyle='-', linewidth=0, markerfacecolor='white', 
                                      markeredgecolor='white', label="// Break"
                                     ),            
                               Line2D([0], [0], color='black',marker='>', markersize=legend_size, 
                                      linestyle='-', linewidth=0, markerfacecolor='white', 
                                      markeredgecolor='black', label="Forward"
                                     ),
                               Line2D([0], [0], color='black',marker='<', markersize=legend_size, 
                                      linestyle='-', linewidth=0, markerfacecolor='white', 
                                      markeredgecolor='black', label="Reverse"
                                     )
                              ]
        else:
            legend_elements = [Line2D([0], [0], color='black',marker='>', markersize=legend_size, 
                                      linestyle='-', linewidth=0, markerfacecolor='white', 
                                      markeredgecolor='black', label="Forward"
                                     ),
                               Line2D([0], [0], color='black',marker='<', markersize=legend_size, 
                                      linestyle='-', linewidth=0, markerfacecolor='white', 
                                      markeredgecolor='black', label="Reverse"
                                     )
                              ]
        
        if colors == None:
            legend_elements.extend(MTColors_legends['MITOFISH'])
            ncol = 12
        elif isinstance(colors, str):
            legend_elements.extend(MTColors_legends.get(colors.upper(), 'MITOFISH'))
            ncol=12
        elif isinstance(colors, dict):
            legend_elements.extend([Patch(facecolor=colors[i], edgecolor='black', label=i) for i in colors if i!="source"])
            ncol = 12
            
        if add_other_genes == False:
            legend_elements = legend_elements[:-1]
        ax.legend(handles=legend_elements,
                  loc='upper right', 
                  bbox_to_anchor=legend_postion, 
                  ncol=12,
                  #frameon=False,
                  #shadow=False,
                  prop={'size': legend_size},
                  title='')
                  
    if output!=None:
        fig.savefig(output, bbox_inches='tight', dpi=dpi)
    if axes != None:
        return ax
    else:
        return fig, ax

def add_tag(axs=None, 
            by_row=True, 
            fontsize=20, 
            fontstyle='bold',
            fontpostiton=(-0.15, 1.1), 
            labels = ["A", "B", "C", "D", "E", "F",
                      "G", "H", "I", "J", "K", "L",
                      "M", "N", "O", "P", "Q", "I",
                      "S", "T", "U", "V", "W", "X",
                      "Y", "Z"]):
    """
    Descripton:
        Add tags to multiple axes.
        
    Parameters：
       axs: {array} a include axes array.
       by_row: {bool}
       fontsize {float}
       fontstyle: {str} font style. such as plain, bold, italic, bold.italic
       fontpostiton: {tuple} font postition.
       labels: a set labels, such as, ["A", "B", "C", ...], ["I", "II", "III", ...].
    
    """
    fontstyles = {"plain": ("normal", "normal"),
     "bold": ("bold", "normal"),
     "italic": ("normal", "italic"),
     "bold.italic": ("bold", "italic")
    }
    
    axs_new= []
    if len(axs.shape) == 1:
        axs_new = axs
    else:
        if by_row:
            for x in range(0, axs.shape[0]):
                for y in range(0, axs.shape[1]):
                    axs_new.append(ax[(x,y)])
        else:
            for x in range(0, axs.shape[1]):
                for y in range(0, axs.shape[0]):
                    axs_new.append(ax[(y,x)])
    for n, ax in enumerate(axs_new):
        
        
        ax.text(*fontpostiton, labels[n], transform=ax.transAxes, fontsize=fontsize, 
                fontweight=fontstyles.get(fontstyle, ("bold", "normal"))[0],
                fontstyle=fontstyles.get(fontstyle, ("bold", "normal"))[1],
                verticalalignment='top')
    return ax

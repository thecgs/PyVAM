# PyVAM
A Python package for visualizing animal mitogenomes.

## Installation

```pip install git+https://github.com/thecgs/PyVAM.git```

## Usage

### Example 1:

You can easily draw a circular mitogenome map in human, just like this, using the [OGDRAW](https://chlorobox.mpimp-golm.mpg.de/OGDraw.html) theme.

```python
import pyvam

fig, ax = pyvam.draw_circos_MT("NC_012920.1", colors="OGDRAW", output="./doc/Fig.1.png", dpi=72)
```

The result is shown in the figure below:

![](doc/Fig.1.png#pic_center)

<center>Fig.1 A circular mitogenome map in human</center>

### Parameter Details:

```python
help(pyvam.draw_circos_MT)

Help on function draw_circos_MT in module pyvam.drawMT:

draw_circos_MT(file, output=None, abbr=False, isfilename2species=False, colors='mitofish', radius=25, show_gene_label=True, gene_label_size=7, gene_label_inner=False, show_info=True, info_fontsize=15, show_legend=True, legend_size=6, legend_postion=(1, -0.1), show_GC_circos=True, GC_circos_height=0.3, GC_circos_color='grey', GC_circos_bin=50, GC_circos_step=50, start=None, axes=None, direction=-1, figsize=(10, 10), tidyname=False, add_id=False, dpi=300, remove_NCR=False, default_topology='circular')
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
                     12S rRNA, 16S rRNA, D-loop.
        axes: {matplotlib.projections.polar.PolarAxes}
        figsize: {tuple} fig of size.
        direction: {1, -1} clockwise: -1, anticlockwise: 1
        tidyname: {bool} tidy gene name.
        add_id: {bool} Species add to accession id from NCBI.
        dpi: {int} dpi value. the resolution in dots per inch.
        remove_NCR: {bool} do not display D-loop.
        default_topology: {str} topology to assume only when the input GenBank
            record does not declare one. Accepts "circular" (default) or
            "linear"; an assumption is reported as a warning.
```



### Example 2:

For mitochondrial genomes that are not closed circles (e.g., due to incomplete assembly), the circular genome map displays a gap at the upper left.

```python
import pyvam

fig, ax = pyvam.draw_circos_MT("MK804157", colors="OGDRAW", output="./doc/Fig.2.png", dpi=72)
```

![](doc/Fig.2.png#pic_center)

<center>Fig.2 A circular mitogenome map in Compsulyx cochereaui</center>

### Example 3:

According to a study by [Wang et al.](https://link.springer.com/article/10.1186/s40850-025-00239-x), compared with the ancestor of Stylommatophora, the mitochondrial genes of *M. pictum* exhibited multiple rearrangement events, while the mitochondrial genes of *S. arundinetorum* showed only minor differences. You can draw the two species separately for comparison (Fig. 3 A–B), or you can draw them together on the same figure, and using the [MitoFish](https://mitofish.aori.u-tokyo.ac.jp/annotation/draw) theme.

```python
import pyvam
import  matplotlib.pyplot as plt

fig, axs = plt.subplots(1, 3, figsize=(20, 20/3), subplot_kw={'projection':'polar'})
plt.subplots_adjust(wspace=0.3)

pyvam.draw_circos_MT(file="OP311642", colors="MitoFish", radius=20, gene_label_size=5, gene_label_inner=False, show_info=True, show_legend=False, axes=axs[0], info_fontsize=6)
pyvam.draw_circos_MT(file="OP289102", colors="MitoFish", radius=20, gene_label_size=5, gene_label_inner=False, show_info=True, show_legend=False, axes=axs[1], info_fontsize=6)
pyvam.draw_circos_MT(file="OP311642", colors="MitoFish", radius=20, gene_label_size=5, gene_label_inner=False, show_info=False,show_legend=False, show_GC_circos=False, axes=axs[2])
pyvam.draw_circos_MT(file="OP289102", colors="MitoFish", radius=12, gene_label_size=4, show_gene_label=True,   show_info=False,show_legend=False, show_GC_circos=False, axes=axs[2])

pyvam.add_tag(axs=axs, by_row=True)

axs[2].text(0.5, 0.5, s="Outer circle: Meghimatium pictum\nInner circle: Succinea arundinetorum", size=6, ha='center', va='center', style='italic')

fig.savefig("./doc/Fig.3.png", bbox_inches='tight', dpi=72)
```

The result is shown in the figure below:

![](doc/Fig.3.png#pic_center)

<center>Fig.3 Three circular mitogenome map in two Stylommatophora species</center>

### Example 4:

Comparing the mitochondrial genomes of multiple species, Circos plots are clearly not the best choice; Below, we use a linear plot to visualize gene rearrangements, and open the circular mitochondrial genome at the specified site (such as, ND1).  [Source of demo data](https://academic.oup.com/isd/article/3/6/12/5686061)

```python
import pyvam

pyvam.draw_linear_MT_nonproportional(files=["MK804148", "MK804158","MK804149", "MK804157"], 
                      start='ND1', add_id=True, dpi=72, force_reoriented=True,
                      output="./doc/Fig.4.png")
```

The result is shown in the figure below:

![](doc/Fig.4.png#pic_center)

<center>Fig.4 Linear plot of nonproportional</center>

### Example 5：

As shown in Fig.4, the genes are not drawn to scale. Below is a scaled drawing of the mitochondrial structure.

```python
import pyvam

pyvam.draw_linear_MT(files=["MK804148", "MK804158","MK804149", "MK804157"], 
                      start='ND1', add_id=True, dpi=72, force_reoriented=True,
                      output="./doc/Fig.5.png")
```

The result is shown in the figure below:

![](doc/Fig.5.png#pic_center)

<center>Fig.5 Linear plot of proportional</center>

### Example 6：

When analyzing large-scale mitochondrial data—including gene overlaps and the identification of intergenic regions—it is necessary to dynamically display the locations of genes. Therefore, pyVAM can also generate dynamic visualizations using Plotly.

For nonproportional:

```python
import pyvam

pyvam.draw_linear_MT_nonproportional_interactive(files=["MK804148", "MK804158", "MK804149", "MK804157"],
                                                  start='COX1', add_id=True, force_reoriented=True,output="./doc/Fig.6.html")
```

[Veiw Fig.6.html](https://html-preview.github.io/?url=https://github.com/thecgs/PyVAM/blob/main/doc/Fig.6.html)

```python
import pyvam

pyvam.draw_linear_MT_interactive(files=["MK804148", "MK804158","MK804149", "MK804157"], 
                                  start='COX1', add_id=True, force_reoriented=True,
                                  output="./doc/Fig.7.html")
```

[Veiw Fig.7.html](https://html-preview.github.io/?url=https://github.com/thecgs/PyVAM/blob/main/doc/Fig.7.html)

## Themes

### Built-in themes:

It comes with 10 built-in themes, including Chen, Tan, OGDRAW, MitoFish, MitoFish1, MitoZ, Chloroplot, Grey, IGV, gggenes.

```python
import pyvam
import matplotlib.pyplot as plt

fig, axs = plt.subplots(2, 5, figsize=(15, 5), subplot_kw={'projection':'polar'})
axes = [a for ax in axs for a in ax]
themes = ["Chen", "Tan", "OGDRAW", "MitoFish", "MitoFish1", "MitoZ", "Chloroplot", "Grey", "IGV", "gggenes"]

for theme, ax in zip(themes, axes):
    pyvam.draw_circos_MT("NC_012920.1", colors=theme, radius=12, show_gene_label=False, 
                              show_info=False, show_legend=False, show_GC_circos=False, axes=ax)
    ax.text(0.5, 0.5, s=theme, ha='center', va='center')
    
fig.savefig("./doc/Fig.8.png", bbox_inches='tight', dpi=72)
```

![](doc/Fig.8.png#pic_center)

<center>Fig.8 10 built-in themes</center>

### Custom colors:

You can also customize the color scheme to your liking. like this,

```python
MTColors_by_Set3 = {'source':"#000000",
                    'ND1': '#8DD3C7',
                    'ND2': '#8DD3C7',
                    'ND3': '#8DD3C7',
                    'ND4L':'#8DD3C7',
                    'ND4': '#8DD3C7',
                    'ND5': '#8DD3C7',
                    'ND6': '#8DD3C7',
                    'COX1':'#BEBADA',
                    'COX2':'#BEBADA',
                    'COX3':'#BEBADA',
                    'ATPase6':'#FDB462',
                    'ATPase8':'#FDB462',
                    'Cytb':'#80B1D3', 
                    'tRNA-His':'#FFFFB3',
                    'tRNA-Pro':'#FFFFB3',
                    'tRNA-Thr':'#FFFFB3',
                    'tRNA-Trp':'#FFFFB3',
                    'tRNA-Met':'#FFFFB3',
                    'tRNA-Asp':'#FFFFB3',
                    'tRNA-Ala':'#FFFFB3',
                    'tRNA-Gln':'#FFFFB3',
                    'tRNA-Ile':'#FFFFB3',
                    'tRNA-Arg':'#FFFFB3',
                    'tRNA-Tyr':'#FFFFB3',
                    'tRNA-Phe':'#FFFFB3',
                    'tRNA-Lys':'#FFFFB3',
                    'tRNA-Gly':'#FFFFB3',
                    'tRNA-Asn':'#FFFFB3',
                    'tRNA-Leu':'#FFFFB3',
                    'tRNA-Glu':'#FFFFB3',
                    'tRNA-Val':'#FFFFB3',
                    'tRNA-Cys':'#FFFFB3',
                    'tRNA-Ser':'#FFFFB3',
                    '12S rRNA':'#FB8072',
                    '16S rRNA':'#FB8072',
                    'D-loop': '#E5C494',
                    'Other genes': '#B3DE69',
                   }
```

## Tidy GenBank

If you're not particularly fond of the PyVAM visualization but still want to use PyVAM's ability to specify the starting point of a GenBank file, you can use this function, which handles the conversion of feature coordinates and the rotation of the mtgenome.

```python
# For the MZ387761, the default starting point is the D-loop.

import pyvam
pyvam.tidy_genbank("MZ387761", 
                    start=None,  # No mitochondrial gene is specified as the starting point.
                    table=2,
                    #output="new_genbank_file.gb"
                   )
```

<details>
    <summary>Genbank file of D-loop starting point </summary>
    <pre><code>
LOCUS       Homo_sapiens                16581 bp    DNA     circular     PRI 08-APR-2026
DEFINITION  .
ACCESSION   .
VERSION     .
KEYWORDS    .
SOURCE      mitochondrion Homo sapiens
  ORGANISM  Homo sapiens
            Unclassified.
REFERENCE   1  (bases 1 to 16581)
  AUTHORS   Chen, G.
  TITLE     PyVAMR: A Python package for visualizing 
            animal mitochondrial rearrangements.
  JOURNAL   Unpublished
  TITLE     Direct Submission
FEATURES             Location/Qualifiers
     source          1..16581
                     /organism="Homo sapiens"
                     /organelle="mitochondrion"
                     /mol_type="genomic DNA"
     D-loop          complement(join(1..580,16036..16581))
                     /note="Control Region"
     tRNA            581..651
                     /product="tRNA-Phe"
     rRNA            652..1605
                     /product="12S ribosomal RNA"
     tRNA            1606..1674
                     /product="tRNA-Val"
     rRNA            1675..3232
                     /product="16S ribosomal RNA"
     tRNA            3233..3307
                     /product="tRNA-Leu"
     gene            3310..4265
                     /gene="ND1"
     CDS             3310..4265
                     /gene="ND1"
                     /codon_start=1
                     /transl_table=2
                     /product="NADH dehydrogenase subunit 1"
                     /translation="MPMANLLLLIVPILIAMAFLMLTERKILGYMQLRKGPNVVGPYG
                     LLQPFADAMKLFTKEPLKPATSTITLYITAPTLALTIALLLWTPLPMPNPLVNLNLGL
                     LFILATSSLAVYSILWSGWASNSNYALIGALRAVAQTISYEVTLAIILLSTLLMSGSF
                     NLSTLITTQEHLWLLLPSWPLAMMWFISTLAETNRTPFDLAEGESELVSGFNIEYAAG
                     PFALFFMAEYTNIIMMNTLTTTIFLGTTYDALSPELYTTYFVTKTLLLTSLFLWIRTA
                     YPRFRYDQLMHLLWKNFLPLTLALLMWHVSMPITISSIPPQT"
     tRNA            4266..4334
                     /product="tRNA-Ile"
     tRNA            complement(4332..4403)
                     /product="tRNA-Gln"
     tRNA            4405..4472
                     /product="tRNA-Met"
     gene            4473..5514
                     /gene="ND2"
     CDS             4473..5514
                     /gene="ND2"
                     /codon_start=1
                     /transl_table=2
                     /product="NADH dehydrogenase subunit 2"
                     /translation="MNPLAQPVIYSTIFAGTLITALSSHWFFTWVGLEMNMLAFIPVL
                     TKKMNPRSTEAAIKYFLTQATASMILLMAILFNNMLSGQWTMTNTTNQYSSLMIMMAM
                     AMKLGMAPFHFWVPEVTQGTPLTSGLLLLTWQKLAPISIMYQISPSLDVSLLLTLSIL
                     SIMAGSWGGLNQTQLRKILAYSSITHMGWMMAVLPYNPNMTILNLTIYIILTTTAFLL
                     LNLNSSTTTLLLSRTWNKLTWLTPLIPSTLLSLGGLPPLTGFLPKWAIIEEFTKNNSL
                     IIPTIMATITLLNLYFYLRLIYSTSITLLPMSNNVKMKWQFEHTKPTPFLPTLIALTT
                     LLLPISPFMLMIL"
     tRNA            5515..5579
                     /product="tRNA-Trp"
     tRNA            complement(5590..5658)
                     /product="tRNA-Ala"
     tRNA            complement(5660..5732)
                     /product="tRNA-Asn"
     tRNA            complement(5764..5829)
                     /product="tRNA-Cys"
     tRNA            complement(5829..5894)
                     /product="tRNA-Tyr"
     gene            5907..7448
                     /gene="COX1"
     CDS             5907..7448
                     /gene="COX1"
                     /codon_start=1
                     /transl_table=2
                     /product="cytochrome c oxidase subunit 1"
                     /translation="MFADRWLFSTNHKDIGTLYLLFGAWAGVLGTALSLLIRAELGQP
                     GNLLGNDHIYNVIVTAHAFVMIFFMVMPIMIGGFGNWLVPLMIGAPDMAFPRMNNMSF
                     WLLPPSLLLLLASAMVETGAGTGWTVYPPLAGNYSHPGASVDLTIFSLHLAGVSSILG
                     AINFITTIINMKPPAMTQYQTPLFVWSVLITAVLLLLSLPVLAAGITMLLTDRNLNTT
                     FFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIVTYYSGKKEPFGYMGMVWA
                     MMSIGFLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAIPTGVKVFSWLATLHGSNMKW
                     SAAVLWALGFIFLFTVGGLTGIVLANSSLDIVLHDTYYVVAHFHYVLSMGAVFAIMGG
                     FIHWFPLFSGYTLDQTYAKIHFTIMFIGVNLTFFPQHFLGLSGMPRRYSDYPDAYTTW
                     NILSSVGSFISLTAVMLMIFMIWEAFASKRKVLMVEEPSMNLEWLYGCPPPYHTFEEP
                     VYMKS"
     tRNA            complement(7448..7519)
                     /product="tRNA-Ser"
     tRNA            7521..7588
                     /product="tRNA-Asp"
     gene            7589..8272
                     /gene="COX2"
     CDS             7589..8272
                     /gene="COX2"
                     /codon_start=1
                     /transl_table=2
                     /product="cytochrome c oxidase subunit 2"
                     /translation="MAHAAQVGLQDATSPIMEELITFHDHALMIIFLICFLVLYALFL
                     TLTTKLTNTNISDAQEMETVWTILPAIILVLIALPSLRILYMTDEVNDPSLTIKSIGH
                     QWYWTYEYTDYGGLIFNSYMLPPLFLEPGDLRLLDVDNRVVLPIEAPIRMMITSQDVL
                     HSWAVPTLGLKTDAIPGRLNQTTFTATRPGVYYGQCSEICGANHSFMPIVLELIPLKI
                     FEMGPVFTL"
     tRNA            8307..8376
                     /product="tRNA-Lys"
     gene            8378..8584
                     /gene="ATPase8"
     CDS             8378..8584
                     /gene="ATPase8"
                     /codon_start=1
                     /transl_table=2
                     /product="ATP synthase F0 subunit 8"
                     /translation="MPQLNTTVWPTMITPMLLTLFLITQLKMLNTNYHLPPSPKPMKM
                     KNYNKPWEPKWTKICSLHSLPPQS"
     gene            8539..9219
                     /gene="ATPase6"
     CDS             8539..9219
                     /gene="ATPase6"
                     /codon_start=1
                     /transl_table=2
                     /product="ATP synthase F0 subunit 6"
                     /translation="MNENLFASFIAPTILGLPAAVLIILFPPLLIPTSKYLINNRLIT
                     TQQWLIKLTSKQMMTMHNTKGRTWSLMLVSLIIFIATTNLLGLLPHSFTPTTQLSMNL
                     AMAIPLWAGAVIMGFRSKIKNALAHFLPQGTPTPLIPMLVIIETISLLIQPMALAVRL
                     TANITAGHLLMHLIGSATLAMSTINLPSTLIIFTILILLTILEIAVALIQAYVFTLLV
                     SLYLHDNT"
     gene            9219..10002
                     /gene="COX3"
     CDS             9219..10002
                     /gene="COX3"
                     /codon_start=1
                     /transl_table=2
                     /product="cytochrome c oxidase subunit 3"
                     /translation="MTHQSHAYHMVKPSPWPLTGALSALLMTSGLAMWFHFHSMTLLM
                     LGLLTNTLTMYQWWRDVTRESTYQGHHTPPVQKGLRYGMILFITSEVFFFAGFFWAFY
                     HSSLAPTPQLGGHWPPTGITPLNPLEVPLLNTSVLLASGVSITWAHHSLMENNRNQMI
                     QALLITILLGLYFTLLQASEYFESPFTISDGIYGSTFFVATGFHGLHVIIGSTFLTIC
                     FIRQLMFHFTSKHHFGFEAAAWYWHFVDVVWLFLYVSIYWWGS"
     tRNA            10003..10070
                     /product="tRNA-Gly"
     gene            10071..10416
                     /gene="ND3"
     CDS             10071..10416
                     /gene="ND3"
                     /codon_start=1
                     /transl_table=2
                     /product="NADH dehydrogenase subunit 3"
                     /translation="MNFALILMINTLLALLLMIITFWLPQLNGYMEKSTPYECGFDPM
                     SPARVPFSMKFFLVAITFLLFDLEIALLLPLPWALQTTNLPLMVMSSLLLIIILALSL
                     AYEWLQKGLDWTE"
     tRNA            10417..10481
                     /product="tRNA-Arg"
     gene            10482..10778
                     /gene="ND4L"
     CDS             10482..10778
                     /gene="ND4L"
                     /codon_start=1
                     /transl_table=2
                     /product="NADH dehydrogenase subunit 4L"
                     /translation="MPLIYMNIMLAFTISLLGMLVYRSHLMSSLLCLEGMMLSLFIMA
                     TLMTLNTHSLLANIVPIAMLVFAACEAAVGLALLVSISNTYGLDYVHNLNLLQC"
     gene            10772..12149
                     /gene="ND4"
     CDS             10772..12149
                     /gene="ND4"
                     /codon_start=1
                     /transl_table=2
                     /product="NADH dehydrogenase subunit 4"
                     /translation="MLKLIVPTIMLLPLTWLSKKHMIWINTTTHSLIISIIPLLFFNQ
                     INNNLFSCSPTFSSDPLTTPLLMLTTWLLPLTIMASQRHLSSEPLSRKKLYLSMLISL
                     QISLIMTFTATELIMFYIFFETTLIPTLAIITRWGNQPERLNAGTYFLFYTLVGSLPL
                     LIALIYTHNTLGSLNILLLTLTAQELSNSWANNLMWLAYTMAFMVKMPLYGLHLWLPK
                     AHVEAPIAGSMVLAAVLLKLGGYGMMRLTLILNPLTKHMAYPFLVLSLWGMIMTSSIC
                     LRQTDLKSLIAYSSISHMALVVTAILIQTPWSFTGAVILMIAHGLTSSLLFCLANSNY
                     ERTHSRIMILSQGLQTLLPLMAFWWLLASLANLALPPTINLLGELSVLVTTFSWSNIT
                     LLLTGLNMLVTALYSLYMFTTTQWGSLTHHINNMKPSFTRENTLMFMHLSPILLLSLN
                     PDIITGFSS"
     tRNA            12150..12218
                     /product="tRNA-His"
     tRNA            12219..12277
                     /product="tRNA-Ser"
     tRNA            12278..12348
                     /product="tRNA-Leu"
     gene            12349..14160
                     /gene="ND5"
     CDS             12349..14160
                     /gene="ND5"
                     /codon_start=1
                     /transl_table=2
                     /product="NADH dehydrogenase subunit 5"
                     /translation="MTMHTTMTTLTLTSLIPPILTTLVNPNKKNSYPHYVKSIVASTF
                     IISLFPTTMFMCLDQEVIISNWHWATTQTTQLSLSFKLDYFSMMFIPVALFVTWSIME
                     FSLWYMNSDPNINQFFKYLLIFLITMLILVTANNLFQLFIGWEGVGIMSFLLISWWYA
                     RADANTAAIQAILYNRIGDIGFILALAWFILHSNSWDPQQMALLNANPSLTPLLGLLL
                     AAAGKSAQLGLHPWLPSAMEGPTPVSALLHSSTMVVAGIFLLIRFHPLAENSPLIQTL
                     TLCLGAITTLFAAVCALTQNDIKKIVAFSTSSQLGLMMVTIGINQPHLAFLHICTHAF
                     FKAMLFMCSGSIIHNLNNEQDIRKMGGLLKTMPLTSTSLTIGSLALAGMPFLTGFYSK
                     DHIIETANMSYTNAWALSITLIATSLTSAYSTRMILLTLTGQPRFPTLTNINENNPTL
                     LNPIKRLAAGSLFAGFLITNNISPASPFQTTIPLYLKLTALAVTFLGLLTALDLNYLT
                     NKLKMKSPLCTFYFSNMLGFYPSITHRTIPYLGLLTSLNLPLLLLDLTWLEKLLPKTI
                     SQHQISTSIITSTQKGMIKLYFLSFFFPLILTLLLIT"
     gene            complement(14161..14685)
                     /gene="ND6"
     CDS             complement(14161..14685)
                     /gene="ND6"
                     /codon_start=1
                     /transl_table=2
                     /product="NADH dehydrogenase subunit 6"
                     /translation="MMYALFLLSVGLVMGFVGFSSKPSPIYGGLVLIVSGVVGCVIIL
                     NFGGGYMGLMVFLIYLGGMMVVFGYTTAMAIEEYPEAWGSGVEVLVSVLVGLAMEVGL
                     VLWVKEYDGVVVVVNFNSVGSWMIYEGEGSGLIREDPIGAGALYDYGRWLVVVTGWTL
                     FVGVYIVIEIARGN"
     tRNA            complement(14686..14754)
                     /product="tRNA-Glu"
     gene            14759..15899
                     /gene="Cytb"
     CDS             14759..15899
                     /gene="Cytb"
                     /codon_start=1
                     /transl_table=2
                     /product="cytochrome b"
                     /translation="MTPMRKINPLMKLINHSFIDLPTPSNISAWWNFGSLLGACLILQ
                     ITTGLFLAMHYSPDASTAFSSIAHITRDVNYGWIIRYLHANGASMFFICLFLHIGRGL
                     YYGSFLYSETWNIGIILLLATMATAFMGYVLPWGQMSFWGATVITNLLSAIPYIGTDL
                     VQWIWGGYSVDSPTLTRFFTFHFILPFIIAALAALHLLFLHETGSNNPLGITSHSDKI
                     TFHPYYTIKDALGLLLFILSLMTLTLFSPDLLGDPDNYTLANPLNTPPHIKPEWYFLF
                     AYTILRSVPNKLGGVLALLLSILILAMIPILHMSKQQSMMFRPLSQSLYWLLAADLLI
                     LTWIGGQPVSYPFTIIGQVASVLYFTTILILMPTISLIENKMLKWA"
     tRNA            15900..15965
                     /product="tRNA-Thr"
     tRNA            complement(15967..16035)
                     /product="tRNA-Pro"
ORIGIN
        1 gatcacaggt ctatcaccct attaaccact cacgggagct ctccatgcat ttggtatttt
       61 cgtctggggg gtgtgcacgc gatagcattg cgagacgctg gagccggagc accctatgtc
      121 gcagtatctg tctttgattc ctgcctcatc ctattattta tcgcacctac gttcaatatt
      181 acaggcgaac atacttacta aagtgtgtta attaattaat gcttgtagga cataataata
      241 acaattgaat gtctgtacag ccgctttcca cacagacatc ataacaaaaa atttccacca
      301 aaccccccct ccccccgctt ctggccacag cacttaaaca catctctgcc aaaccccaaa
      361 aacaaagaac cctaacacca gcctaaccag atttcaaatt ttatcttttg gcggtatgca
      421 cttttaacag tcacccccca actaacacat tattttcccc tcccactccc atactactaa
      481 tctcatcaat acaacccccg cccatcctac ccagcacaca cacaccgctg ctaaccccat
      541 accccgaacc aaccaaaccc caaagacacc cccccccaca gtttatgtag cttacctcct
      601 caaagcaata cactgaaaat gtttagacgg gctcacatca ccccataaac aaataggttt
      661 ggtcctagcc tttctattag ctcttagtaa gattacacat gcaagcatcc ccattccagt
      721 gagttcaccc tctaaatcac cacgatcaaa agggacaagc atcaagcacg cagcaatgca
      781 gctcaaaacg cttagcctag ccacaccccc acgggaaaca gcagtgatta acctttagca
      841 ataaacgaaa gtttaactaa gctatactaa ccccagggtt ggtcaatttc gtgccagcca
      901 ccgcggtcac acgattaacc caagtcaata gaagccggcg taaagagtgt tttagatcac
      961 cccctcccca ataaagctaa aactcacctg agttgtaaaa aactccagtt gacacaaaat
     1021 agactacgaa agtggcttta acatatctga acacacaata gctaagaccc aaactgggat
     1081 tagatacccc actatgctta gccctaaacc tcaacagtta aatcaacaaa actgctcgcc
     1141 agaacactac gagccacagc ttaaaactca aaggacctgg cggtgcttca tatccctcta
     1201 gaggagcctg ttctgtaatc gataaacccc gatcaacctc accacctctt gctcagccta
     1261 tataccgcca tcttcagcaa accctgatga aggctacaaa gtaagcgcaa gtacccacgt
     1321 aaagacgtta ggtcaaggtg tagcccatga ggtggcaaga aatgggctac attttctacc
     1381 ccagaaaact acgatagccc ttatgaaact taagggtcga aggtggattt agcagtaaac
     1441 tgagagtaga gtgcttagtt gaacagggcc ctgaagcgcg tacacaccgc ccgtcaccct
     1501 cctcaagtat acttcaaagg acatttaact aaaaccccta cgcatttata tagaggagac
     1561 aagtcgtaac atggtaagtg tactggaaag tgcacttgga cgaaccagag tgtagcttaa
     1621 cacaaagcac ccaacttaca cttaggagat ttcaacttaa cttgaccgct ctgagctaaa
     1681 cctagcccca aacccactcc accttactac cagacaacct tagccaaacc atttacccaa
     1741 ataaagtata ggcgatagaa attgaaacct ggcgcaatag atatagtacc gcaagggaaa
     1801 gatgaaaaat tataaccaag cataatatag caaggactaa cccctatacc ttctgcataa
     1861 tgaattaact agaaataact ttgcaaggag aaccaaagct aagacccccg aaaccagacg
     1921 agctacctaa gaacagctaa aagagcacac ccgtctatgt agcaaaatag tgggaagatt
     1981 tataggtaga ggcgacaaac ctaccgagcc tggtgatagc tggttgtcca agatagaatc
     2041 ttagttcaac tttaaatttg cccacagaac cctctaaatc cccttgtaaa tttaactgtt
     2101 agtccaaaga ggaacagctc tttggacact aggaaaaaac cttgtagaga gagtaaaaaa
     2161 tttaacaccc atagtaggcc taaaagcagc caccaattaa gaaagcgttc aagctcaaca
     2221 cccactacct aaaaaatccc aaacatataa ctgaactcct cacacccaat tggaccaatc
     2281 tatcacccta tagaagaact aatgttagta taagtaacat gaaaacattc tcctccgcat
     2341 aagcctgcgt cagattaaaa cactgaactg acaattaaca gcccaatatc tacaatcaac
     2401 caacaagtca ttattaccct cactgtcaac ccaacacagg catgctcata aggaaaggtt
     2461 aaaaaaagta aaaggaactc ggcaaatctt accccgcctg tttaccaaaa acatcacctc
     2521 tagcatcacc agtattagag gcaccgcctg cccagtgaca catgtttaac ggccgcggta
     2581 ccctaaccgt gcaaaggtag cataatcact tgttccttaa atagggacct gtatgaatgg
     2641 ctccacgagg gttcagctgt ctcttacttt taaccagtga aattgacctg cccgtgaaga
     2701 ggcgggcatg acacagcaag acgagaagac cctatggagc tttaatttat taatgcaaac
     2761 agtacctaac aaacccacag gtcctaaact accaaacctg cattaaaaat ttcggttggg
     2821 gcgacctcgg agcagaaccc aacctccgag cagtacatgc taagacttca ccagtcaaag
     2881 cgaactacta tactcaattg atccaataac ttgaccaacg gaacaagtta ccctagggat
     2941 aacagcgcaa tcctattcta gagtccatat caacaatagg gtttacgacc tcgatgttgg
     3001 atcaggacat cccgatggtg cagccgctat taaaggttcg tttgttcaac gattaaagtc
     3061 ctacgtgatc tgagttcaga ccggagtaat ccaggtcggt ttctatctac ttcaaattcc
     3121 tccctgtacg aaaggacaag agaaataagg cctacttcac aaagcgcctt cccccgtaaa
     3181 tgatatcatc tcaacttagt attataccca cacccaccca agaacagggt ttgttaagat
     3241 ggcagagccc ggtaatcgca taaaacttaa aactttacag tcagaggttc aattcctctt
     3301 cttaacaaca tacccatggc caacctccta ctcctcattg tacccattct aatcgcaatg
     3361 gcattcctaa tgcttaccga acgaaaaatt ctaggctata tacaactacg caaaggcccc
     3421 aacgttgtag gcccctacgg gctactacaa cccttcgctg acgccataaa actcttcacc
     3481 aaagagcccc taaaacccgc cacatctacc atcaccctct acatcaccgc cccgacctta
     3541 gctctcacca tcgctcttct actatgaacc cccctcccca tacccaaccc cctggtcaac
     3601 ctcaacctag gcctcctatt tattctagcc acctctagcc tagccgttta ctcaatcctc
     3661 tgatcagggt gagcatcaaa ctcaaactac gccctgatcg gcgcactgcg agcagtagcc
     3721 caaacaatct catatgaagt caccctagcc atcattctac tatcaacatt actaataagt
     3781 ggctccttta acctctccac ccttatcaca acacaagaac acctctgatt actcctgcca
     3841 tcatgaccct tggccataat atgatttatc tccacactag cagagaccaa ccgaaccccc
     3901 ttcgaccttg ccgaagggga gtccgaacta gtctcaggct tcaacatcga atacgccgca
     3961 ggccccttcg ccctattctt catagccgaa tacacaaata ttattataat aaacaccctc
     4021 accactacaa tcttcctagg aacaacatat gacgcactct cccctgaact ctacacaaca
     4081 tattttgtca ccaagaccct acttctaacc tccctgttct tatgaattcg aacagcatac
     4141 ccccgattcc gctacgacca actcatacac ctcctatgaa aaaacttcct accactcacc
     4201 ctagcattac ttatatgaca tgtctccata cccattacaa tctccagcat tccccctcaa
     4261 acctaagaaa tatgtctgat aaaagagtta ctttgataga gtaaataata ggagcttaaa
     4321 cccccttatt tctaggacta tgagaatcga acccatccct gagaatccaa aattctccgt
     4381 gccacctatc acaccccatc ctaaagtaag gtcagctaaa taagctatcg ggcccatacc
     4441 ccgaaaatgt tggttatacc cttcccgtac taattaatcc cctggcccaa cccgtcatct
     4501 actctaccat ctttgcaggc acactcatca cagcgctaag ctcgcactga ttttttacct
     4561 gagtaggcct agaaataaac atgctagctt ttattccagt tctaaccaaa aaaataaacc
     4621 ctcgttccac agaagctgcc atcaagtatt tcctcacgca agcaaccgca tccataatcc
     4681 ttctaatagc tatcctcttc aacaatatac tctccggaca atgaaccata accaatacta
     4741 ccaatcaata ctcatcatta ataatcataa tggctatagc aataaaacta ggaatagccc
     4801 cctttcactt ctgagtccca gaggttaccc aaggcacccc tctgacatcc ggcctgcttc
     4861 ttctcacatg acaaaaacta gcccccatct caatcatata ccaaatctct ccctcactag
     4921 acgtaagcct tctcctcact ctctcaatct tatccatcat agcaggcagt tgaggtggat
     4981 taaaccaaac ccagctacgc aaaatcttag catactcctc aattacccac ataggatgaa
     5041 taatagcagt tctaccgtac aaccctaaca taaccattct taatttaact atttatatta
     5101 tcctaactac taccgcattc ctactactca acttaaactc cagcaccacg accctactac
     5161 tatctcgcac ctgaaacaag ctaacatgac taacaccctt aattccatcc accctcctct
     5221 ccctaggagg cctgcccccg ctaaccggct ttttgcccaa atgggccatt atcgaagaat
     5281 tcacaaaaaa caatagcctc atcatcccca ccatcatagc caccatcacc ctccttaacc
     5341 tctacttcta cctacgccta atctactcca cctcaatcac actactcccc atatctaaca
     5401 acgtaaaaat aaaatgacag tttgaacata caaaacccac cccattcctc cccacactca
     5461 tcgcccttac cacgctactc ctacctatct ccccttttat actaataatc ttatagaaat
     5521 ttaggttaaa tacagaccaa gagccttcaa agccctcagt aagttgcaat acttaatttc
     5581 tgtaacagct aaggactgca aaaccccact ctgcatcaac tgaacgcaaa tcagccactt
     5641 taattaagct aagcccttac tagaccaatg ggacttaaac ccacaaacac ttagttaaca
     5701 gctaagcacc ctaatcaact ggcttcaatc tacttctccc gccgccggga aaaaaggcgg
     5761 gagaagcccc ggcaggtttg aagctgcttc ttcgaatttg caattcaata tgaaaatcat
     5821 ctcggagctg gtaaaaagag gcctaacccc tgtctttaga tttacagtcc aatgcttcac
     5881 tcagccattt tacctcaccc ccactgatgt tcgccgaccg ttgactattc tctacaaacc
     5941 acaaagacat tggaacacta tacctattat tcggcgcatg agctggagtc ctaggcacag
     6001 ctctaagcct ccttattcga gccgagctgg gccagccagg caaccttcta ggtaacgacc
     6061 acatctacaa cgttatcgtc acagcccatg catttgtaat aatcttcttc atagtaatac
     6121 ccatcataat cggaggcttt ggcaactgac tagttcccct aataatcggt gcccccgata
     6181 tggcgtttcc ccgcataaac aacataagct tctgactctt acctccctct ctcctactcc
     6241 tgctcgcatc tgctatagtg gagaccggag caggaacagg ttgaacagtc taccctccct
     6301 tagcagggaa ctactcccac cctggagcct ccgtagacct aaccatcttc tccttacacc
     6361 tagcaggtgt ctcctctatc ttaggggcca tcaatttcat cacaacaatt atcaatataa
     6421 aaccccctgc cataacccaa taccaaacgc ccctcttcgt ctgatccgtc ctaatcacag
     6481 cagtcctact tctcctatct ctcccagtcc tagctgctgg catcactata ctactaacag
     6541 accgcaacct caacaccacc ttcttcgacc ccgccggagg aggagacccc attctatacc
     6601 aacacctatt ctgatttttc ggtcaccctg aagtttatat tcttatccta ccaggcttcg
     6661 gaataatctc ccatattgta acttactact ccggaaaaaa agaaccattt ggatacatag
     6721 gtatggtctg agctatgata tcaattggct tcctagggtt tatcgtgtga gcacaccata
     6781 tatttacagt aggaatagac gtagacacac gagcatattt cacctccgct accataatca
     6841 tcgctatccc caccggcgtc aaagtattta gctgactcgc cacactccac ggaagcaata
     6901 tgaaatgatc tgctgcagtg ctctgagccc taggattcat ctttcttttc accgtaggtg
     6961 gcctgactgg cattgtatta gcaaactcat cactagacat cgtactacac gacacgtact
     7021 acgttgtagc tcacttccac tatgtcctat caataggagc tgtatttgcc atcataggag
     7081 gcttcattca ctgatttccc ctattctcag gctacaccct agaccaaacc tacgccaaaa
     7141 tccatttcac tatcatattc atcggcgtaa atctaacttt cttcccacaa cactttctcg
     7201 gcctatccgg aatgccccga cgttactcgg actaccccga tgcatacacc acatgaaaca
     7261 tcctatcatc tgtaggctca ttcatttctc taacagcagt aatattaata attttcatga
     7321 tttgagaagc cttcgcttcg aagcgaaaag tcctaatagt agaagaaccc tccataaacc
     7381 tggagtgact atatggatgc cccccaccct accacacatt cgaagaaccc gtatacataa
     7441 aatctagaca aaaaaggaag gaatcgaacc ccccaaagct ggtttcaagc caaccccatg
     7501 gcctccatga ctttttcaaa aaggtattag aaaaaccatt tcataacttt gtcaaagtta
     7561 aattataggc taaatcctat atatcttaat ggcacatgca gcgcaagtag gtctacaaga
     7621 cgctacttcc cctatcatag aagagcttat cacctttcat gatcacgccc tcataatcat
     7681 tttccttatc tgcttcctag tcctgtatgc ccttttccta acactcacaa caaaactaac
     7741 taatactaac atctcagacg ctcaggaaat agaaaccgtc tgaactatcc tgcccgccat
     7801 catcctagtc ctcatcgccc tcccatccct acgcatcctt tacataacag acgaggtcaa
     7861 cgatccctcc cttaccatca aatcaattgg ccaccaatgg tactgaacct acgagtacac
     7921 cgactacggc ggactaatct tcaactccta catacttccc ccattattcc tagaaccagg
     7981 cgacctgcga ctccttgacg ttgacaatcg agtagtactc ccgattgaag cccccattcg
     8041 tataataatt acatcacaag acgtcttgca ctcatgagct gtccccacat taggcttaaa
     8101 aacagatgca attcccggac gtctaaacca aaccactttc accgctacac gaccgggggt
     8161 atactacggt caatgctctg aaatctgtgg agcaaaccac agtttcatgc ccatcgtcct
     8221 agaattaatt cccctaaaaa tctttgaaat agggcccgta tttaccctat agcaccccct
     8281 ctaccccctc taccccctct agagcccact gtaaagctaa cttagcatta accttttaag
     8341 ttaaagatta agagaaccaa cacctcttta cagtgaaatg ccccaactaa atactaccgt
     8401 atggcccacc ataattaccc ccatactcct tacactattc ctcatcaccc aactaaaaat
     8461 attaaataca aactaccacc tacctccctc accaaagccc ataaaaataa aaaattataa
     8521 caaaccctga gaaccaaaat gaacgaaaat ctgttcgctt cattcattgc ccccacaatc
     8581 ctaggcctac ccgccgcagt actgatcatt ctatttcccc ctctattgat ccccacctcc
     8641 aaatatctca tcaacaaccg actaatcacc acccaacaat gactaatcaa actaacctca
     8701 aaacaaataa taaccataca caacactaaa ggacgaacct gatctcttat actagtatcc
     8761 ttaatcattt ttattgccac aactaacctc ctcggactcc tgcctcactc atttacacca
     8821 accacccaac tatctataaa cctagccatg gccatcccct tatgagcggg cgcagtgatt
     8881 ataggctttc gctctaagat taaaaatgcc ctagcccact tcttaccaca aggcacacct
     8941 acacccctta tccccatact agttattatc gaaaccatca gcctactcat tcaaccaata
     9001 gccctggccg tacgcctaac cgctaacatt actgcaggcc acctactcat gcacctaatt
     9061 ggaagcgcca ccctagcaat atcaaccatt aaccttccct ctacacttat catcttcaca
     9121 attctaattc tactgactat cctagaaatc gctgtcgcct taatccaagc ctacgttttc
     9181 acacttctag taagcctcta cctgcacgac aacacataat gacccaccaa tcacatgcct
     9241 atcatatagt aaaacccagc ccatgacccc taacaggggc cctctcagcc ctcctaatga
     9301 cctccggcct agccatgtga tttcacttcc actccataac gctcctcata ctaggcctac
     9361 taaccaacac actaaccata taccaatgat ggcgcgatgt aacacgagaa agcacatacc
     9421 aaggccacca cacaccacct gtccaaaaag gccttcgata cgggataatc ctatttatta
     9481 cctcagaagt ttttttcttc gcaggatttt tctgagcctt ttaccactcc agcctagccc
     9541 ctacccccca attaggaggg cactggcccc caacaggcat caccccgcta aatcccctag
     9601 aagtcccact cctaaacaca tccgtattac tcgcatcagg agtatcaatc acctgagctc
     9661 accatagtct aatagaaaac aaccgaaacc aaataattca agcactgctt attacaattt
     9721 tactgggtct ctattttacc ctcctacaag cctcagagta cttcgagtct cccttcacca
     9781 tttccgacgg catctacggc tcaacatttt ttgtagccac aggcttccac ggacttcacg
     9841 tcattattgg ctcaactttc ctcactatct gcttcatccg ccaactaata tttcacttta
     9901 catccaaaca tcactttggc ttcgaagccg ccgcctgata ctggcatttt gtagatgtgg
     9961 tttgactatt tctgtatgtc tccatctatt gatgagggtc ttactctttt agtataaata
    10021 gtaccgttaa cttccaatta actagttttg acaacattca aaaaagagta ataaacttcg
    10081 ccttaatttt aataatcaac accctcctag ccttactact aataattatt acattttgac
    10141 taccacaact caacggctac atagaaaaat ccacccctta cgagtgcggc ttcgacccta
    10201 tatcccccgc ccgcgtccct ttctccataa aattcttctt agtagctatt accttcttat
    10261 tatttgatct agaaattgcc ctccttttac ccctaccatg agccctacaa acaactaacc
    10321 tgccactaat agttatgtca tccctcttat taatcatcat cctagcccta agtctggcct
    10381 atgagtgact acaaaaagga ttagactgaa ccgaattggt atatagttta aacaaaacga
    10441 atgatttcga ctcattaaat tatgataatc atatctacca aatgcccctc atttacataa
    10501 atattatact agcatttacc atctcacttc taggaatact agtatatcgc tcacacctca
    10561 tatcctccct actatgccta gaaggaataa tactatcgct gttcattata gctactctca
    10621 taaccctcaa cacccactcc ctcttagcca atattgtgcc tattgccata ctagtctttg
    10681 ccgcctgcga agcagcggtg ggcctagccc tactagtctc aatctccaac acatatggcc
    10741 tagactacgt acataaccta aacctactcc aatgctaaaa ctaatcgtcc caacaattat
    10801 attactacca ctgacatgac tttccaaaaa acatataatt tgaatcaaca caaccaccca
    10861 cagcctaatt attagcatca tccctctact attttttaac caaatcaaca acaacctatt
    10921 tagctgttcc ccaacctttt cctccgaccc cctaacaacc cccctcctaa tactaactac
    10981 ctgactccta cccctcacaa tcatggcaag ccaacgccac ttatccagtg aaccactatc
    11041 acgaaaaaaa ctctacctct ctatactaat ctccctacaa atctccttaa ttataacatt
    11101 cacagccaca gaactaatca tattttatat cttcttcgaa accacactta tccccacctt
    11161 ggctatcatc acccgatgag gcaaccagcc agaacgcctg aacgcaggca catacttcct
    11221 attctacacc ctagtaggct cccttcccct actcatcgca ctgatttaca ctcacaacac
    11281 cctaggctca ctaaacattc tactactcac tctcactgcc caagaactat caaactcctg
    11341 agccaacaac ttaatatgac tagcttacac aatagctttt atagtaaaga tacctcttta
    11401 cggactccac ttatgactcc ctaaagccca tgtcgaagcc cccatcgctg ggtcaatagt
    11461 acttgccgca gtactcttaa aactaggcgg ctatggtata atacgcctca cactcattct
    11521 caaccccctg acaaaacaca tagcctaccc cttccttgta ctatccctat gaggcataat
    11581 tataacaagc tccatctgcc tacgacaaac agacctaaaa tcgctcattg catactcttc
    11641 aatcagccac atagccctcg tagtaacagc cattctcatc caaaccccct gaagcttcac
    11701 cggcgcagtc attctcataa tcgcccacgg acttacatcc tcattactat tctgcctagc
    11761 aaactcaaac tacgaacgca ctcacagtcg catcataatc ctctctcaag gacttcaaac
    11821 tctgctccca ctaatagctt tttgatgact tctagcaagc ctcgctaacc tcgccttacc
    11881 ccccactatt aacctactgg gagaactctc tgtgctagta accacgttct cctgatcaaa
    11941 tatcactctc ctacttacag gactcaacat actagtcaca gccctatact ccctctacat
    12001 atttaccaca acacaatggg gctcactcac ccaccacatt aacaacataa aaccctcatt
    12061 cacacgagaa aacaccctca tgttcataca cctatccccc attctcctcc tatccctcaa
    12121 ccccgacatc attaccgggt tttcctcttg taaatatagt ttaaccaaaa catcagattg
    12181 tgaatctgac aacagaggct tacgacccct tatttaccga gaaagctcac aagaactgct
    12241 aactcatgcc cccatgtcta acaacatggc tttctcaact tttaaaggat aacagctatc
    12301 cattggtctt aggccccaaa aattttggtg caactccaaa taaaagtaat aaccatgcac
    12361 actactataa ccaccctaac cctgacttcc ctaattcccc ccatccttac caccctcgtt
    12421 aaccctaaca aaaaaaactc atacccccat tatgtaaaat ccattgtcgc atccaccttt
    12481 attatcagtc tcttccccac aacaatattc atgtgcctag accaagaagt tattatctcg
    12541 aactgacact gagccacaac ccaaacaacc cagctctccc taagcttcaa actagactac
    12601 ttctccataa tattcatccc tgtagcattg ttcgttacat ggtccatcat agaattctca
    12661 ctgtgatata taaactcaga cccaaacatt aatcagttct tcaaatatct actcatcttc
    12721 ctaattacca tactaatctt agttaccgct aacaacctat tccaactgtt catcggctga
    12781 gagggcgtag gaattatatc cttcttgctc atcagttgat gatacgcccg agcagatgcc
    12841 aacacagcag ccattcaagc aatcctatac aaccgtatcg gcgatatcgg tttcatcctc
    12901 gccttagcat gatttatcct acactccaac tcatgagacc cacaacaaat agcccttcta
    12961 aacgctaatc caagcctcac cccactacta ggcctcctcc tagcagcagc aggcaaatca
    13021 gcccaattag gtctccaccc ctgactcccc tcagccatag aaggccccac cccagtctca
    13081 gccctactcc actcaagcac tatagttgta gcaggaatct tcttactcat ccgcttccac
    13141 cccctagcag aaaatagccc actaatccaa actctaacac tatgcttagg cgctatcacc
    13201 actctgttcg cagcagtctg cgcccttaca caaaatgaca tcaaaaagat cgtagccttc
    13261 tccacttcaa gtcaactagg actcataata gttacaatcg gcatcaacca accacaccta
    13321 gcattcctgc acatctgtac ccacgccttc ttcaaagcca tactatttat gtgctccgga
    13381 tccatcatcc acaaccttaa caatgaacaa gatattcgaa aaataggagg actactcaaa
    13441 accatacctc tcacttcaac ctccctcacc attggcagcc tagcattagc aggaatacct
    13501 ttcctcacag gtttctactc caaagaccac atcatcgaaa ccgcaaacat atcatacaca
    13561 aacgcctgag ccctatctat tactctcatc gctacctccc tgacaagcgc ctatagcact
    13621 cgaataattc ttctcaccct aacaggtcaa cctcgcttcc ccacccttac taacattaac
    13681 gaaaataacc ccaccctact aaaccccatt aaacgcctgg cagccggaag cctattcgca
    13741 ggatttctca ttactaacaa catttccccc gcatccccct tccaaacaac aatccccctc
    13801 tacctaaaac tcacagccct cgctgtcact ttcctaggac ttctaacagc cctagacctc
    13861 aactacctaa ccaacaaact taaaataaaa tccccactat gcacatttta tttctccaac
    13921 atactcggat tctaccctag catcacacac cgcacaatcc cctatctagg ccttcttacg
    13981 agcctaaacc tacccctact cctcctagac ctaacctgac tagaaaagct attacctaaa
    14041 acaatttcac agcaccaaat ctccacctcc atcatcacct caacccaaaa aggcataatt
    14101 aaactttact tcctctcttt cttcttccca ctcatcctaa ccctactcct aatcacataa
    14161 cctattcccc cgagcaatct caattacaat atatacacca acaaacaatg ttcaaccagt
    14221 aactactact aatcaacgcc catagtcata caaagccccc gcaccaatag gatcctcccg
    14281 aatcaaccct gacccctctc cttcataaat tattcagctt cctacactat taaagtttac
    14341 cacaaccacc accccatcat actctttcac ccacagcacc aatcctacct ccatcgctaa
    14401 ccccactaaa acactcacca agacctcaac ccctgacccc catgcctcag gatactcctc
    14461 aatagccatc gctgtagtat atccaaagac aaccatcatt ccccctaaat aaattaaaaa
    14521 aactattaaa cccatataac ctcccccaaa attcagaata ataacacacc cgaccacacc
    14581 gctaacaatc aatactaaac ccccataaat aggagaaggc ttagaagaaa accccacaaa
    14641 ccccattact aaacccacac tcaacagaaa caaagcatac atcattattc tcgcacggac
    14701 tacaaccacg accaatgata tgaaaaacca tcgttgtatt tcaactacaa gaacaccaat
    14761 gaccccaata cgcaaaatta accccctaat aaaattaatt aaccactcat tcatcgacct
    14821 ccccacccca tccaacatct ccgcatgatg aaacttcggc tcactccttg gcgcctgcct
    14881 gatcctccaa atcaccacag gactattcct agccatacac tactcaccag acgcctcaac
    14941 cgccttttca tcaatcgccc acatcactcg agacgtaaat tatggctgaa tcatccgcta
    15001 ccttcacgcc aatggcgcct caatattctt tatctgcctc ttcctacaca tcgggcgagg
    15061 cctatattac ggatcatttc tctactcaga aacctgaaac atcggcatta tcctcctgct
    15121 tgcaactata gcaacagcct tcataggcta tgtcctcccg tgaggccaaa tatcattctg
    15181 aggggccaca gtaattacaa acttactatc cgccatccca tacattggga cagacctagt
    15241 tcaatgaatc tgaggaggct actcagtaga cagtcccacc ctcacacgat tctttacctt
    15301 tcacttcatc ttgcccttca ttattgcagc cctagcagca ctccacctcc tattcttgca
    15361 cgaaacggga tcaaacaacc ccctaggaat cacctcccat tccgataaaa tcaccttcca
    15421 cccttactac acaatcaaag acgccctcgg cttacttctc ttcattctct ccttaatgac
    15481 attaacacta ttctcaccag acctcctagg cgacccagac aattataccc tagccaaccc
    15541 cttaaacacc cctccccaca tcaagcccga atgatatttc ctattcgcct acacaattct
    15601 ccgatccgtc cctaacaagc taggaggcgt ccttgcccta ttactatcca tcctcatcct
    15661 agcaataatc cccatcctcc atatatccaa acaacaaagc ataatatttc gcccactaag
    15721 ccaatcactt tattgactcc tagccgcaga cctcctcatt ctaacctgaa tcggaggaca
    15781 accagtaagc taccctttta ccatcattgg acaagtagca tccgtactat acttcacaac
    15841 aatcctaatc ctaataccaa ctatctccct aattgaaaac aaaatactca aatgggcctg
    15901 tccttgtagt ataaactaat acaccagtct tgtaaaccga agatgaaaac ctttttccaa
    15961 ggacaaatca gagaaaaagt ctttaactcc accattagca cccaaagcta agattctaat
    16021 ttaaactatt ctctgttctt tcatggggaa gcagatttgg gtaccaccca agtattgact
    16081 cacccatcaa caaccgctat gtatttcgta cattactgcc agccaccatg aatattgcac
    16141 ggtaccataa atacttgacc acctgtagta cataaaaacc caatccacat caaaaccccc
    16201 tccccatgct tacaagcaag tacagcaatc aaccctcaac tatcacacat caactgcaac
    16261 tccaaagcca cccctcaccc actaggatac caacaaacct acctatcctt aacagtacat
    16321 agtacataaa gccatttacc gtacatagca cattacagtc aaatcccttc tcgtccccat
    16381 ggatgacccc cctcagatag gggtcccttg accaccatcc tccgtgaaat caatatcccg
    16441 cacaagagtg ctactctcct cgctccgggc ccataacact tgggggtagc taaagtgaac
    16501 tgtatccgac atctggttcc tacttcaggg ccataaagcc taaatagccc acacgttccc
    16561 cttaaataag acatcacgat g
//
</code></pre> 
</details>


```python
# For MZ387761, specify tRNA-Phe as the starting point.

import pyvam
pyvam.tidy_genbank("MZ387761", 
                    start='tRNA-Phe',
                    table=2,
                    #output="new_genbank_file.gb"
                   )
```

<details>
    <summary> Genbank of tRNA-Phe starting point </summary>
    <pre><code>
LOCUS       Homo_sapiens                16581 bp    DNA     circular     PRI 08-APR-2026
DEFINITION  .
ACCESSION   .
VERSION     .
KEYWORDS    .
SOURCE      mitochondrion Homo sapiens
  ORGANISM  Homo sapiens
            Unclassified.
REFERENCE   1  (bases 1 to 16581)
  AUTHORS   Chen, G.
  TITLE     PyVAMR: A Python package for visualizing 
            animal mitochondrial rearrangements.
  JOURNAL   Unpublished
  TITLE     Direct Submission
FEATURES             Location/Qualifiers
     source          1..16581
                     /organism="Homo sapiens"
                     /organelle="mitochondrion"
                     /mol_type="genomic DNA"
     tRNA            1..71
                     /product="tRNA-Phe"
     rRNA            72..1025
                     /product="12S ribosomal RNA"
     tRNA            1026..1094
                     /product="tRNA-Val"
     rRNA            1095..2652
                     /product="16S ribosomal RNA"
     tRNA            2653..2727
                     /product="tRNA-Leu"
     gene            2730..3685
                     /gene="ND1"
     CDS             2730..3685
                     /gene="ND1"
                     /codon_start=1
                     /transl_table=2
                     /product="NADH dehydrogenase subunit 1"
                     /translation="MPMANLLLLIVPILIAMAFLMLTERKILGYMQLRKGPNVVGPYG
                     LLQPFADAMKLFTKEPLKPATSTITLYITAPTLALTIALLLWTPLPMPNPLVNLNLGL
                     LFILATSSLAVYSILWSGWASNSNYALIGALRAVAQTISYEVTLAIILLSTLLMSGSF
                     NLSTLITTQEHLWLLLPSWPLAMMWFISTLAETNRTPFDLAEGESELVSGFNIEYAAG
                     PFALFFMAEYTNIIMMNTLTTTIFLGTTYDALSPELYTTYFVTKTLLLTSLFLWIRTA
                     YPRFRYDQLMHLLWKNFLPLTLALLMWHVSMPITISSIPPQT"
     tRNA            3686..3754
                     /product="tRNA-Ile"
     tRNA            complement(3752..3823)
                     /product="tRNA-Gln"
     tRNA            3825..3892
                     /product="tRNA-Met"
     gene            3893..4934
                     /gene="ND2"
     CDS             3893..4934
                     /gene="ND2"
                     /codon_start=1
                     /transl_table=2
                     /product="NADH dehydrogenase subunit 2"
                     /translation="MNPLAQPVIYSTIFAGTLITALSSHWFFTWVGLEMNMLAFIPVL
                     TKKMNPRSTEAAIKYFLTQATASMILLMAILFNNMLSGQWTMTNTTNQYSSLMIMMAM
                     AMKLGMAPFHFWVPEVTQGTPLTSGLLLLTWQKLAPISIMYQISPSLDVSLLLTLSIL
                     SIMAGSWGGLNQTQLRKILAYSSITHMGWMMAVLPYNPNMTILNLTIYIILTTTAFLL
                     LNLNSSTTTLLLSRTWNKLTWLTPLIPSTLLSLGGLPPLTGFLPKWAIIEEFTKNNSL
                     IIPTIMATITLLNLYFYLRLIYSTSITLLPMSNNVKMKWQFEHTKPTPFLPTLIALTT
                     LLLPISPFMLMIL"
     tRNA            4935..4999
                     /product="tRNA-Trp"
     tRNA            complement(5010..5078)
                     /product="tRNA-Ala"
     tRNA            complement(5080..5152)
                     /product="tRNA-Asn"
     tRNA            complement(5184..5249)
                     /product="tRNA-Cys"
     tRNA            complement(5249..5314)
                     /product="tRNA-Tyr"
     gene            5327..6868
                     /gene="COX1"
     CDS             5327..6868
                     /gene="COX1"
                     /codon_start=1
                     /transl_table=2
                     /product="cytochrome c oxidase subunit 1"
                     /translation="MFADRWLFSTNHKDIGTLYLLFGAWAGVLGTALSLLIRAELGQP
                     GNLLGNDHIYNVIVTAHAFVMIFFMVMPIMIGGFGNWLVPLMIGAPDMAFPRMNNMSF
                     WLLPPSLLLLLASAMVETGAGTGWTVYPPLAGNYSHPGASVDLTIFSLHLAGVSSILG
                     AINFITTIINMKPPAMTQYQTPLFVWSVLITAVLLLLSLPVLAAGITMLLTDRNLNTT
                     FFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIVTYYSGKKEPFGYMGMVWA
                     MMSIGFLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAIPTGVKVFSWLATLHGSNMKW
                     SAAVLWALGFIFLFTVGGLTGIVLANSSLDIVLHDTYYVVAHFHYVLSMGAVFAIMGG
                     FIHWFPLFSGYTLDQTYAKIHFTIMFIGVNLTFFPQHFLGLSGMPRRYSDYPDAYTTW
                     NILSSVGSFISLTAVMLMIFMIWEAFASKRKVLMVEEPSMNLEWLYGCPPPYHTFEEP
                     VYMKS"
     tRNA            complement(6868..6939)
                     /product="tRNA-Ser"
     tRNA            6941..7008
                     /product="tRNA-Asp"
     gene            7009..7692
                     /gene="COX2"
     CDS             7009..7692
                     /gene="COX2"
                     /codon_start=1
                     /transl_table=2
                     /product="cytochrome c oxidase subunit 2"
                     /translation="MAHAAQVGLQDATSPIMEELITFHDHALMIIFLICFLVLYALFL
                     TLTTKLTNTNISDAQEMETVWTILPAIILVLIALPSLRILYMTDEVNDPSLTIKSIGH
                     QWYWTYEYTDYGGLIFNSYMLPPLFLEPGDLRLLDVDNRVVLPIEAPIRMMITSQDVL
                     HSWAVPTLGLKTDAIPGRLNQTTFTATRPGVYYGQCSEICGANHSFMPIVLELIPLKI
                     FEMGPVFTL"
     tRNA            7727..7796
                     /product="tRNA-Lys"
     gene            7798..8004
                     /gene="ATPase8"
     CDS             7798..8004
                     /gene="ATPase8"
                     /codon_start=1
                     /transl_table=2
                     /product="ATP synthase F0 subunit 8"
                     /translation="MPQLNTTVWPTMITPMLLTLFLITQLKMLNTNYHLPPSPKPMKM
                     KNYNKPWEPKWTKICSLHSLPPQS"
     gene            7959..8639
                     /gene="ATPase6"
     CDS             7959..8639
                     /gene="ATPase6"
                     /codon_start=1
                     /transl_table=2
                     /product="ATP synthase F0 subunit 6"
                     /translation="MNENLFASFIAPTILGLPAAVLIILFPPLLIPTSKYLINNRLIT
                     TQQWLIKLTSKQMMTMHNTKGRTWSLMLVSLIIFIATTNLLGLLPHSFTPTTQLSMNL
                     AMAIPLWAGAVIMGFRSKIKNALAHFLPQGTPTPLIPMLVIIETISLLIQPMALAVRL
                     TANITAGHLLMHLIGSATLAMSTINLPSTLIIFTILILLTILEIAVALIQAYVFTLLV
                     SLYLHDNT"
     gene            8639..9422
                     /gene="COX3"
     CDS             8639..9422
                     /gene="COX3"
                     /codon_start=1
                     /transl_table=2
                     /product="cytochrome c oxidase subunit 3"
                     /translation="MTHQSHAYHMVKPSPWPLTGALSALLMTSGLAMWFHFHSMTLLM
                     LGLLTNTLTMYQWWRDVTRESTYQGHHTPPVQKGLRYGMILFITSEVFFFAGFFWAFY
                     HSSLAPTPQLGGHWPPTGITPLNPLEVPLLNTSVLLASGVSITWAHHSLMENNRNQMI
                     QALLITILLGLYFTLLQASEYFESPFTISDGIYGSTFFVATGFHGLHVIIGSTFLTIC
                     FIRQLMFHFTSKHHFGFEAAAWYWHFVDVVWLFLYVSIYWWGS"
     tRNA            9423..9490
                     /product="tRNA-Gly"
     gene            9491..9836
                     /gene="ND3"
     CDS             9491..9836
                     /gene="ND3"
                     /codon_start=1
                     /transl_table=2
                     /product="NADH dehydrogenase subunit 3"
                     /translation="MNFALILMINTLLALLLMIITFWLPQLNGYMEKSTPYECGFDPM
                     SPARVPFSMKFFLVAITFLLFDLEIALLLPLPWALQTTNLPLMVMSSLLLIIILALSL
                     AYEWLQKGLDWTE"
     tRNA            9837..9901
                     /product="tRNA-Arg"
     gene            9902..10198
                     /gene="ND4L"
     CDS             9902..10198
                     /gene="ND4L"
                     /codon_start=1
                     /transl_table=2
                     /product="NADH dehydrogenase subunit 4L"
                     /translation="MPLIYMNIMLAFTISLLGMLVYRSHLMSSLLCLEGMMLSLFIMA
                     TLMTLNTHSLLANIVPIAMLVFAACEAAVGLALLVSISNTYGLDYVHNLNLLQC"
     gene            10192..11569
                     /gene="ND4"
     CDS             10192..11569
                     /gene="ND4"
                     /codon_start=1
                     /transl_table=2
                     /product="NADH dehydrogenase subunit 4"
                     /translation="MLKLIVPTIMLLPLTWLSKKHMIWINTTTHSLIISIIPLLFFNQ
                     INNNLFSCSPTFSSDPLTTPLLMLTTWLLPLTIMASQRHLSSEPLSRKKLYLSMLISL
                     QISLIMTFTATELIMFYIFFETTLIPTLAIITRWGNQPERLNAGTYFLFYTLVGSLPL
                     LIALIYTHNTLGSLNILLLTLTAQELSNSWANNLMWLAYTMAFMVKMPLYGLHLWLPK
                     AHVEAPIAGSMVLAAVLLKLGGYGMMRLTLILNPLTKHMAYPFLVLSLWGMIMTSSIC
                     LRQTDLKSLIAYSSISHMALVVTAILIQTPWSFTGAVILMIAHGLTSSLLFCLANSNY
                     ERTHSRIMILSQGLQTLLPLMAFWWLLASLANLALPPTINLLGELSVLVTTFSWSNIT
                     LLLTGLNMLVTALYSLYMFTTTQWGSLTHHINNMKPSFTRENTLMFMHLSPILLLSLN
                     PDIITGFSS"
     tRNA            11570..11638
                     /product="tRNA-His"
     tRNA            11639..11697
                     /product="tRNA-Ser"
     tRNA            11698..11768
                     /product="tRNA-Leu"
     gene            11769..13580
                     /gene="ND5"
     CDS             11769..13580
                     /gene="ND5"
                     /codon_start=1
                     /transl_table=2
                     /product="NADH dehydrogenase subunit 5"
                     /translation="MTMHTTMTTLTLTSLIPPILTTLVNPNKKNSYPHYVKSIVASTF
                     IISLFPTTMFMCLDQEVIISNWHWATTQTTQLSLSFKLDYFSMMFIPVALFVTWSIME
                     FSLWYMNSDPNINQFFKYLLIFLITMLILVTANNLFQLFIGWEGVGIMSFLLISWWYA
                     RADANTAAIQAILYNRIGDIGFILALAWFILHSNSWDPQQMALLNANPSLTPLLGLLL
                     AAAGKSAQLGLHPWLPSAMEGPTPVSALLHSSTMVVAGIFLLIRFHPLAENSPLIQTL
                     TLCLGAITTLFAAVCALTQNDIKKIVAFSTSSQLGLMMVTIGINQPHLAFLHICTHAF
                     FKAMLFMCSGSIIHNLNNEQDIRKMGGLLKTMPLTSTSLTIGSLALAGMPFLTGFYSK
                     DHIIETANMSYTNAWALSITLIATSLTSAYSTRMILLTLTGQPRFPTLTNINENNPTL
                     LNPIKRLAAGSLFAGFLITNNISPASPFQTTIPLYLKLTALAVTFLGLLTALDLNYLT
                     NKLKMKSPLCTFYFSNMLGFYPSITHRTIPYLGLLTSLNLPLLLLDLTWLEKLLPKTI
                     SQHQISTSIITSTQKGMIKLYFLSFFFPLILTLLLIT"
     gene            complement(13581..14105)
                     /gene="ND6"
     CDS             complement(13581..14105)
                     /gene="ND6"
                     /codon_start=1
                     /transl_table=2
                     /product="NADH dehydrogenase subunit 6"
                     /translation="MMYALFLLSVGLVMGFVGFSSKPSPIYGGLVLIVSGVVGCVIIL
                     NFGGGYMGLMVFLIYLGGMMVVFGYTTAMAIEEYPEAWGSGVEVLVSVLVGLAMEVGL
                     VLWVKEYDGVVVVVNFNSVGSWMIYEGEGSGLIREDPIGAGALYDYGRWLVVVTGWTL
                     FVGVYIVIEIARGN"
     tRNA            complement(14106..14174)
                     /product="tRNA-Glu"
     gene            14179..15319
                     /gene="Cytb"
     CDS             14179..15319
                     /gene="Cytb"
                     /codon_start=1
                     /transl_table=2
                     /product="cytochrome b"
                     /translation="MTPMRKINPLMKLINHSFIDLPTPSNISAWWNFGSLLGACLILQ
                     ITTGLFLAMHYSPDASTAFSSIAHITRDVNYGWIIRYLHANGASMFFICLFLHIGRGL
                     YYGSFLYSETWNIGIILLLATMATAFMGYVLPWGQMSFWGATVITNLLSAIPYIGTDL
                     VQWIWGGYSVDSPTLTRFFTFHFILPFIIAALAALHLLFLHETGSNNPLGITSHSDKI
                     TFHPYYTIKDALGLLLFILSLMTLTLFSPDLLGDPDNYTLANPLNTPPHIKPEWYFLF
                     AYTILRSVPNKLGGVLALLLSILILAMIPILHMSKQQSMMFRPLSQSLYWLLAADLLI
                     LTWIGGQPVSYPFTIIGQVASVLYFTTILILMPTISLIENKMLKWA"
     tRNA            15320..15385
                     /product="tRNA-Thr"
     tRNA            complement(15387..15455)
                     /product="tRNA-Pro"
     D-loop          complement(15456..16581)
                     /note="Control Region"
ORIGIN
        1 gtttatgtag cttacctcct caaagcaata cactgaaaat gtttagacgg gctcacatca
       61 ccccataaac aaataggttt ggtcctagcc tttctattag ctcttagtaa gattacacat
      121 gcaagcatcc ccattccagt gagttcaccc tctaaatcac cacgatcaaa agggacaagc
      181 atcaagcacg cagcaatgca gctcaaaacg cttagcctag ccacaccccc acgggaaaca
      241 gcagtgatta acctttagca ataaacgaaa gtttaactaa gctatactaa ccccagggtt
      301 ggtcaatttc gtgccagcca ccgcggtcac acgattaacc caagtcaata gaagccggcg
      361 taaagagtgt tttagatcac cccctcccca ataaagctaa aactcacctg agttgtaaaa
      421 aactccagtt gacacaaaat agactacgaa agtggcttta acatatctga acacacaata
      481 gctaagaccc aaactgggat tagatacccc actatgctta gccctaaacc tcaacagtta
      541 aatcaacaaa actgctcgcc agaacactac gagccacagc ttaaaactca aaggacctgg
      601 cggtgcttca tatccctcta gaggagcctg ttctgtaatc gataaacccc gatcaacctc
      661 accacctctt gctcagccta tataccgcca tcttcagcaa accctgatga aggctacaaa
      721 gtaagcgcaa gtacccacgt aaagacgtta ggtcaaggtg tagcccatga ggtggcaaga
      781 aatgggctac attttctacc ccagaaaact acgatagccc ttatgaaact taagggtcga
      841 aggtggattt agcagtaaac tgagagtaga gtgcttagtt gaacagggcc ctgaagcgcg
      901 tacacaccgc ccgtcaccct cctcaagtat acttcaaagg acatttaact aaaaccccta
      961 cgcatttata tagaggagac aagtcgtaac atggtaagtg tactggaaag tgcacttgga
     1021 cgaaccagag tgtagcttaa cacaaagcac ccaacttaca cttaggagat ttcaacttaa
     1081 cttgaccgct ctgagctaaa cctagcccca aacccactcc accttactac cagacaacct
     1141 tagccaaacc atttacccaa ataaagtata ggcgatagaa attgaaacct ggcgcaatag
     1201 atatagtacc gcaagggaaa gatgaaaaat tataaccaag cataatatag caaggactaa
     1261 cccctatacc ttctgcataa tgaattaact agaaataact ttgcaaggag aaccaaagct
     1321 aagacccccg aaaccagacg agctacctaa gaacagctaa aagagcacac ccgtctatgt
     1381 agcaaaatag tgggaagatt tataggtaga ggcgacaaac ctaccgagcc tggtgatagc
     1441 tggttgtcca agatagaatc ttagttcaac tttaaatttg cccacagaac cctctaaatc
     1501 cccttgtaaa tttaactgtt agtccaaaga ggaacagctc tttggacact aggaaaaaac
     1561 cttgtagaga gagtaaaaaa tttaacaccc atagtaggcc taaaagcagc caccaattaa
     1621 gaaagcgttc aagctcaaca cccactacct aaaaaatccc aaacatataa ctgaactcct
     1681 cacacccaat tggaccaatc tatcacccta tagaagaact aatgttagta taagtaacat
     1741 gaaaacattc tcctccgcat aagcctgcgt cagattaaaa cactgaactg acaattaaca
     1801 gcccaatatc tacaatcaac caacaagtca ttattaccct cactgtcaac ccaacacagg
     1861 catgctcata aggaaaggtt aaaaaaagta aaaggaactc ggcaaatctt accccgcctg
     1921 tttaccaaaa acatcacctc tagcatcacc agtattagag gcaccgcctg cccagtgaca
     1981 catgtttaac ggccgcggta ccctaaccgt gcaaaggtag cataatcact tgttccttaa
     2041 atagggacct gtatgaatgg ctccacgagg gttcagctgt ctcttacttt taaccagtga
     2101 aattgacctg cccgtgaaga ggcgggcatg acacagcaag acgagaagac cctatggagc
     2161 tttaatttat taatgcaaac agtacctaac aaacccacag gtcctaaact accaaacctg
     2221 cattaaaaat ttcggttggg gcgacctcgg agcagaaccc aacctccgag cagtacatgc
     2281 taagacttca ccagtcaaag cgaactacta tactcaattg atccaataac ttgaccaacg
     2341 gaacaagtta ccctagggat aacagcgcaa tcctattcta gagtccatat caacaatagg
     2401 gtttacgacc tcgatgttgg atcaggacat cccgatggtg cagccgctat taaaggttcg
     2461 tttgttcaac gattaaagtc ctacgtgatc tgagttcaga ccggagtaat ccaggtcggt
     2521 ttctatctac ttcaaattcc tccctgtacg aaaggacaag agaaataagg cctacttcac
     2581 aaagcgcctt cccccgtaaa tgatatcatc tcaacttagt attataccca cacccaccca
     2641 agaacagggt ttgttaagat ggcagagccc ggtaatcgca taaaacttaa aactttacag
     2701 tcagaggttc aattcctctt cttaacaaca tacccatggc caacctccta ctcctcattg
     2761 tacccattct aatcgcaatg gcattcctaa tgcttaccga acgaaaaatt ctaggctata
     2821 tacaactacg caaaggcccc aacgttgtag gcccctacgg gctactacaa cccttcgctg
     2881 acgccataaa actcttcacc aaagagcccc taaaacccgc cacatctacc atcaccctct
     2941 acatcaccgc cccgacctta gctctcacca tcgctcttct actatgaacc cccctcccca
     3001 tacccaaccc cctggtcaac ctcaacctag gcctcctatt tattctagcc acctctagcc
     3061 tagccgttta ctcaatcctc tgatcagggt gagcatcaaa ctcaaactac gccctgatcg
     3121 gcgcactgcg agcagtagcc caaacaatct catatgaagt caccctagcc atcattctac
     3181 tatcaacatt actaataagt ggctccttta acctctccac ccttatcaca acacaagaac
     3241 acctctgatt actcctgcca tcatgaccct tggccataat atgatttatc tccacactag
     3301 cagagaccaa ccgaaccccc ttcgaccttg ccgaagggga gtccgaacta gtctcaggct
     3361 tcaacatcga atacgccgca ggccccttcg ccctattctt catagccgaa tacacaaata
     3421 ttattataat aaacaccctc accactacaa tcttcctagg aacaacatat gacgcactct
     3481 cccctgaact ctacacaaca tattttgtca ccaagaccct acttctaacc tccctgttct
     3541 tatgaattcg aacagcatac ccccgattcc gctacgacca actcatacac ctcctatgaa
     3601 aaaacttcct accactcacc ctagcattac ttatatgaca tgtctccata cccattacaa
     3661 tctccagcat tccccctcaa acctaagaaa tatgtctgat aaaagagtta ctttgataga
     3721 gtaaataata ggagcttaaa cccccttatt tctaggacta tgagaatcga acccatccct
     3781 gagaatccaa aattctccgt gccacctatc acaccccatc ctaaagtaag gtcagctaaa
     3841 taagctatcg ggcccatacc ccgaaaatgt tggttatacc cttcccgtac taattaatcc
     3901 cctggcccaa cccgtcatct actctaccat ctttgcaggc acactcatca cagcgctaag
     3961 ctcgcactga ttttttacct gagtaggcct agaaataaac atgctagctt ttattccagt
     4021 tctaaccaaa aaaataaacc ctcgttccac agaagctgcc atcaagtatt tcctcacgca
     4081 agcaaccgca tccataatcc ttctaatagc tatcctcttc aacaatatac tctccggaca
     4141 atgaaccata accaatacta ccaatcaata ctcatcatta ataatcataa tggctatagc
     4201 aataaaacta ggaatagccc cctttcactt ctgagtccca gaggttaccc aaggcacccc
     4261 tctgacatcc ggcctgcttc ttctcacatg acaaaaacta gcccccatct caatcatata
     4321 ccaaatctct ccctcactag acgtaagcct tctcctcact ctctcaatct tatccatcat
     4381 agcaggcagt tgaggtggat taaaccaaac ccagctacgc aaaatcttag catactcctc
     4441 aattacccac ataggatgaa taatagcagt tctaccgtac aaccctaaca taaccattct
     4501 taatttaact atttatatta tcctaactac taccgcattc ctactactca acttaaactc
     4561 cagcaccacg accctactac tatctcgcac ctgaaacaag ctaacatgac taacaccctt
     4621 aattccatcc accctcctct ccctaggagg cctgcccccg ctaaccggct ttttgcccaa
     4681 atgggccatt atcgaagaat tcacaaaaaa caatagcctc atcatcccca ccatcatagc
     4741 caccatcacc ctccttaacc tctacttcta cctacgccta atctactcca cctcaatcac
     4801 actactcccc atatctaaca acgtaaaaat aaaatgacag tttgaacata caaaacccac
     4861 cccattcctc cccacactca tcgcccttac cacgctactc ctacctatct ccccttttat
     4921 actaataatc ttatagaaat ttaggttaaa tacagaccaa gagccttcaa agccctcagt
     4981 aagttgcaat acttaatttc tgtaacagct aaggactgca aaaccccact ctgcatcaac
     5041 tgaacgcaaa tcagccactt taattaagct aagcccttac tagaccaatg ggacttaaac
     5101 ccacaaacac ttagttaaca gctaagcacc ctaatcaact ggcttcaatc tacttctccc
     5161 gccgccggga aaaaaggcgg gagaagcccc ggcaggtttg aagctgcttc ttcgaatttg
     5221 caattcaata tgaaaatcat ctcggagctg gtaaaaagag gcctaacccc tgtctttaga
     5281 tttacagtcc aatgcttcac tcagccattt tacctcaccc ccactgatgt tcgccgaccg
     5341 ttgactattc tctacaaacc acaaagacat tggaacacta tacctattat tcggcgcatg
     5401 agctggagtc ctaggcacag ctctaagcct ccttattcga gccgagctgg gccagccagg
     5461 caaccttcta ggtaacgacc acatctacaa cgttatcgtc acagcccatg catttgtaat
     5521 aatcttcttc atagtaatac ccatcataat cggaggcttt ggcaactgac tagttcccct
     5581 aataatcggt gcccccgata tggcgtttcc ccgcataaac aacataagct tctgactctt
     5641 acctccctct ctcctactcc tgctcgcatc tgctatagtg gagaccggag caggaacagg
     5701 ttgaacagtc taccctccct tagcagggaa ctactcccac cctggagcct ccgtagacct
     5761 aaccatcttc tccttacacc tagcaggtgt ctcctctatc ttaggggcca tcaatttcat
     5821 cacaacaatt atcaatataa aaccccctgc cataacccaa taccaaacgc ccctcttcgt
     5881 ctgatccgtc ctaatcacag cagtcctact tctcctatct ctcccagtcc tagctgctgg
     5941 catcactata ctactaacag accgcaacct caacaccacc ttcttcgacc ccgccggagg
     6001 aggagacccc attctatacc aacacctatt ctgatttttc ggtcaccctg aagtttatat
     6061 tcttatccta ccaggcttcg gaataatctc ccatattgta acttactact ccggaaaaaa
     6121 agaaccattt ggatacatag gtatggtctg agctatgata tcaattggct tcctagggtt
     6181 tatcgtgtga gcacaccata tatttacagt aggaatagac gtagacacac gagcatattt
     6241 cacctccgct accataatca tcgctatccc caccggcgtc aaagtattta gctgactcgc
     6301 cacactccac ggaagcaata tgaaatgatc tgctgcagtg ctctgagccc taggattcat
     6361 ctttcttttc accgtaggtg gcctgactgg cattgtatta gcaaactcat cactagacat
     6421 cgtactacac gacacgtact acgttgtagc tcacttccac tatgtcctat caataggagc
     6481 tgtatttgcc atcataggag gcttcattca ctgatttccc ctattctcag gctacaccct
     6541 agaccaaacc tacgccaaaa tccatttcac tatcatattc atcggcgtaa atctaacttt
     6601 cttcccacaa cactttctcg gcctatccgg aatgccccga cgttactcgg actaccccga
     6661 tgcatacacc acatgaaaca tcctatcatc tgtaggctca ttcatttctc taacagcagt
     6721 aatattaata attttcatga tttgagaagc cttcgcttcg aagcgaaaag tcctaatagt
     6781 agaagaaccc tccataaacc tggagtgact atatggatgc cccccaccct accacacatt
     6841 cgaagaaccc gtatacataa aatctagaca aaaaaggaag gaatcgaacc ccccaaagct
     6901 ggtttcaagc caaccccatg gcctccatga ctttttcaaa aaggtattag aaaaaccatt
     6961 tcataacttt gtcaaagtta aattataggc taaatcctat atatcttaat ggcacatgca
     7021 gcgcaagtag gtctacaaga cgctacttcc cctatcatag aagagcttat cacctttcat
     7081 gatcacgccc tcataatcat tttccttatc tgcttcctag tcctgtatgc ccttttccta
     7141 acactcacaa caaaactaac taatactaac atctcagacg ctcaggaaat agaaaccgtc
     7201 tgaactatcc tgcccgccat catcctagtc ctcatcgccc tcccatccct acgcatcctt
     7261 tacataacag acgaggtcaa cgatccctcc cttaccatca aatcaattgg ccaccaatgg
     7321 tactgaacct acgagtacac cgactacggc ggactaatct tcaactccta catacttccc
     7381 ccattattcc tagaaccagg cgacctgcga ctccttgacg ttgacaatcg agtagtactc
     7441 ccgattgaag cccccattcg tataataatt acatcacaag acgtcttgca ctcatgagct
     7501 gtccccacat taggcttaaa aacagatgca attcccggac gtctaaacca aaccactttc
     7561 accgctacac gaccgggggt atactacggt caatgctctg aaatctgtgg agcaaaccac
     7621 agtttcatgc ccatcgtcct agaattaatt cccctaaaaa tctttgaaat agggcccgta
     7681 tttaccctat agcaccccct ctaccccctc taccccctct agagcccact gtaaagctaa
     7741 cttagcatta accttttaag ttaaagatta agagaaccaa cacctcttta cagtgaaatg
     7801 ccccaactaa atactaccgt atggcccacc ataattaccc ccatactcct tacactattc
     7861 ctcatcaccc aactaaaaat attaaataca aactaccacc tacctccctc accaaagccc
     7921 ataaaaataa aaaattataa caaaccctga gaaccaaaat gaacgaaaat ctgttcgctt
     7981 cattcattgc ccccacaatc ctaggcctac ccgccgcagt actgatcatt ctatttcccc
     8041 ctctattgat ccccacctcc aaatatctca tcaacaaccg actaatcacc acccaacaat
     8101 gactaatcaa actaacctca aaacaaataa taaccataca caacactaaa ggacgaacct
     8161 gatctcttat actagtatcc ttaatcattt ttattgccac aactaacctc ctcggactcc
     8221 tgcctcactc atttacacca accacccaac tatctataaa cctagccatg gccatcccct
     8281 tatgagcggg cgcagtgatt ataggctttc gctctaagat taaaaatgcc ctagcccact
     8341 tcttaccaca aggcacacct acacccctta tccccatact agttattatc gaaaccatca
     8401 gcctactcat tcaaccaata gccctggccg tacgcctaac cgctaacatt actgcaggcc
     8461 acctactcat gcacctaatt ggaagcgcca ccctagcaat atcaaccatt aaccttccct
     8521 ctacacttat catcttcaca attctaattc tactgactat cctagaaatc gctgtcgcct
     8581 taatccaagc ctacgttttc acacttctag taagcctcta cctgcacgac aacacataat
     8641 gacccaccaa tcacatgcct atcatatagt aaaacccagc ccatgacccc taacaggggc
     8701 cctctcagcc ctcctaatga cctccggcct agccatgtga tttcacttcc actccataac
     8761 gctcctcata ctaggcctac taaccaacac actaaccata taccaatgat ggcgcgatgt
     8821 aacacgagaa agcacatacc aaggccacca cacaccacct gtccaaaaag gccttcgata
     8881 cgggataatc ctatttatta cctcagaagt ttttttcttc gcaggatttt tctgagcctt
     8941 ttaccactcc agcctagccc ctacccccca attaggaggg cactggcccc caacaggcat
     9001 caccccgcta aatcccctag aagtcccact cctaaacaca tccgtattac tcgcatcagg
     9061 agtatcaatc acctgagctc accatagtct aatagaaaac aaccgaaacc aaataattca
     9121 agcactgctt attacaattt tactgggtct ctattttacc ctcctacaag cctcagagta
     9181 cttcgagtct cccttcacca tttccgacgg catctacggc tcaacatttt ttgtagccac
     9241 aggcttccac ggacttcacg tcattattgg ctcaactttc ctcactatct gcttcatccg
     9301 ccaactaata tttcacttta catccaaaca tcactttggc ttcgaagccg ccgcctgata
     9361 ctggcatttt gtagatgtgg tttgactatt tctgtatgtc tccatctatt gatgagggtc
     9421 ttactctttt agtataaata gtaccgttaa cttccaatta actagttttg acaacattca
     9481 aaaaagagta ataaacttcg ccttaatttt aataatcaac accctcctag ccttactact
     9541 aataattatt acattttgac taccacaact caacggctac atagaaaaat ccacccctta
     9601 cgagtgcggc ttcgacccta tatcccccgc ccgcgtccct ttctccataa aattcttctt
     9661 agtagctatt accttcttat tatttgatct agaaattgcc ctccttttac ccctaccatg
     9721 agccctacaa acaactaacc tgccactaat agttatgtca tccctcttat taatcatcat
     9781 cctagcccta agtctggcct atgagtgact acaaaaagga ttagactgaa ccgaattggt
     9841 atatagttta aacaaaacga atgatttcga ctcattaaat tatgataatc atatctacca
     9901 aatgcccctc atttacataa atattatact agcatttacc atctcacttc taggaatact
     9961 agtatatcgc tcacacctca tatcctccct actatgccta gaaggaataa tactatcgct
    10021 gttcattata gctactctca taaccctcaa cacccactcc ctcttagcca atattgtgcc
    10081 tattgccata ctagtctttg ccgcctgcga agcagcggtg ggcctagccc tactagtctc
    10141 aatctccaac acatatggcc tagactacgt acataaccta aacctactcc aatgctaaaa
    10201 ctaatcgtcc caacaattat attactacca ctgacatgac tttccaaaaa acatataatt
    10261 tgaatcaaca caaccaccca cagcctaatt attagcatca tccctctact attttttaac
    10321 caaatcaaca acaacctatt tagctgttcc ccaacctttt cctccgaccc cctaacaacc
    10381 cccctcctaa tactaactac ctgactccta cccctcacaa tcatggcaag ccaacgccac
    10441 ttatccagtg aaccactatc acgaaaaaaa ctctacctct ctatactaat ctccctacaa
    10501 atctccttaa ttataacatt cacagccaca gaactaatca tattttatat cttcttcgaa
    10561 accacactta tccccacctt ggctatcatc acccgatgag gcaaccagcc agaacgcctg
    10621 aacgcaggca catacttcct attctacacc ctagtaggct cccttcccct actcatcgca
    10681 ctgatttaca ctcacaacac cctaggctca ctaaacattc tactactcac tctcactgcc
    10741 caagaactat caaactcctg agccaacaac ttaatatgac tagcttacac aatagctttt
    10801 atagtaaaga tacctcttta cggactccac ttatgactcc ctaaagccca tgtcgaagcc
    10861 cccatcgctg ggtcaatagt acttgccgca gtactcttaa aactaggcgg ctatggtata
    10921 atacgcctca cactcattct caaccccctg acaaaacaca tagcctaccc cttccttgta
    10981 ctatccctat gaggcataat tataacaagc tccatctgcc tacgacaaac agacctaaaa
    11041 tcgctcattg catactcttc aatcagccac atagccctcg tagtaacagc cattctcatc
    11101 caaaccccct gaagcttcac cggcgcagtc attctcataa tcgcccacgg acttacatcc
    11161 tcattactat tctgcctagc aaactcaaac tacgaacgca ctcacagtcg catcataatc
    11221 ctctctcaag gacttcaaac tctgctccca ctaatagctt tttgatgact tctagcaagc
    11281 ctcgctaacc tcgccttacc ccccactatt aacctactgg gagaactctc tgtgctagta
    11341 accacgttct cctgatcaaa tatcactctc ctacttacag gactcaacat actagtcaca
    11401 gccctatact ccctctacat atttaccaca acacaatggg gctcactcac ccaccacatt
    11461 aacaacataa aaccctcatt cacacgagaa aacaccctca tgttcataca cctatccccc
    11521 attctcctcc tatccctcaa ccccgacatc attaccgggt tttcctcttg taaatatagt
    11581 ttaaccaaaa catcagattg tgaatctgac aacagaggct tacgacccct tatttaccga
    11641 gaaagctcac aagaactgct aactcatgcc cccatgtcta acaacatggc tttctcaact
    11701 tttaaaggat aacagctatc cattggtctt aggccccaaa aattttggtg caactccaaa
    11761 taaaagtaat aaccatgcac actactataa ccaccctaac cctgacttcc ctaattcccc
    11821 ccatccttac caccctcgtt aaccctaaca aaaaaaactc atacccccat tatgtaaaat
    11881 ccattgtcgc atccaccttt attatcagtc tcttccccac aacaatattc atgtgcctag
    11941 accaagaagt tattatctcg aactgacact gagccacaac ccaaacaacc cagctctccc
    12001 taagcttcaa actagactac ttctccataa tattcatccc tgtagcattg ttcgttacat
    12061 ggtccatcat agaattctca ctgtgatata taaactcaga cccaaacatt aatcagttct
    12121 tcaaatatct actcatcttc ctaattacca tactaatctt agttaccgct aacaacctat
    12181 tccaactgtt catcggctga gagggcgtag gaattatatc cttcttgctc atcagttgat
    12241 gatacgcccg agcagatgcc aacacagcag ccattcaagc aatcctatac aaccgtatcg
    12301 gcgatatcgg tttcatcctc gccttagcat gatttatcct acactccaac tcatgagacc
    12361 cacaacaaat agcccttcta aacgctaatc caagcctcac cccactacta ggcctcctcc
    12421 tagcagcagc aggcaaatca gcccaattag gtctccaccc ctgactcccc tcagccatag
    12481 aaggccccac cccagtctca gccctactcc actcaagcac tatagttgta gcaggaatct
    12541 tcttactcat ccgcttccac cccctagcag aaaatagccc actaatccaa actctaacac
    12601 tatgcttagg cgctatcacc actctgttcg cagcagtctg cgcccttaca caaaatgaca
    12661 tcaaaaagat cgtagccttc tccacttcaa gtcaactagg actcataata gttacaatcg
    12721 gcatcaacca accacaccta gcattcctgc acatctgtac ccacgccttc ttcaaagcca
    12781 tactatttat gtgctccgga tccatcatcc acaaccttaa caatgaacaa gatattcgaa
    12841 aaataggagg actactcaaa accatacctc tcacttcaac ctccctcacc attggcagcc
    12901 tagcattagc aggaatacct ttcctcacag gtttctactc caaagaccac atcatcgaaa
    12961 ccgcaaacat atcatacaca aacgcctgag ccctatctat tactctcatc gctacctccc
    13021 tgacaagcgc ctatagcact cgaataattc ttctcaccct aacaggtcaa cctcgcttcc
    13081 ccacccttac taacattaac gaaaataacc ccaccctact aaaccccatt aaacgcctgg
    13141 cagccggaag cctattcgca ggatttctca ttactaacaa catttccccc gcatccccct
    13201 tccaaacaac aatccccctc tacctaaaac tcacagccct cgctgtcact ttcctaggac
    13261 ttctaacagc cctagacctc aactacctaa ccaacaaact taaaataaaa tccccactat
    13321 gcacatttta tttctccaac atactcggat tctaccctag catcacacac cgcacaatcc
    13381 cctatctagg ccttcttacg agcctaaacc tacccctact cctcctagac ctaacctgac
    13441 tagaaaagct attacctaaa acaatttcac agcaccaaat ctccacctcc atcatcacct
    13501 caacccaaaa aggcataatt aaactttact tcctctcttt cttcttccca ctcatcctaa
    13561 ccctactcct aatcacataa cctattcccc cgagcaatct caattacaat atatacacca
    13621 acaaacaatg ttcaaccagt aactactact aatcaacgcc catagtcata caaagccccc
    13681 gcaccaatag gatcctcccg aatcaaccct gacccctctc cttcataaat tattcagctt
    13741 cctacactat taaagtttac cacaaccacc accccatcat actctttcac ccacagcacc
    13801 aatcctacct ccatcgctaa ccccactaaa acactcacca agacctcaac ccctgacccc
    13861 catgcctcag gatactcctc aatagccatc gctgtagtat atccaaagac aaccatcatt
    13921 ccccctaaat aaattaaaaa aactattaaa cccatataac ctcccccaaa attcagaata
    13981 ataacacacc cgaccacacc gctaacaatc aatactaaac ccccataaat aggagaaggc
    14041 ttagaagaaa accccacaaa ccccattact aaacccacac tcaacagaaa caaagcatac
    14101 atcattattc tcgcacggac tacaaccacg accaatgata tgaaaaacca tcgttgtatt
    14161 tcaactacaa gaacaccaat gaccccaata cgcaaaatta accccctaat aaaattaatt
    14221 aaccactcat tcatcgacct ccccacccca tccaacatct ccgcatgatg aaacttcggc
    14281 tcactccttg gcgcctgcct gatcctccaa atcaccacag gactattcct agccatacac
    14341 tactcaccag acgcctcaac cgccttttca tcaatcgccc acatcactcg agacgtaaat
    14401 tatggctgaa tcatccgcta ccttcacgcc aatggcgcct caatattctt tatctgcctc
    14461 ttcctacaca tcgggcgagg cctatattac ggatcatttc tctactcaga aacctgaaac
    14521 atcggcatta tcctcctgct tgcaactata gcaacagcct tcataggcta tgtcctcccg
    14581 tgaggccaaa tatcattctg aggggccaca gtaattacaa acttactatc cgccatccca
    14641 tacattggga cagacctagt tcaatgaatc tgaggaggct actcagtaga cagtcccacc
    14701 ctcacacgat tctttacctt tcacttcatc ttgcccttca ttattgcagc cctagcagca
    14761 ctccacctcc tattcttgca cgaaacggga tcaaacaacc ccctaggaat cacctcccat
    14821 tccgataaaa tcaccttcca cccttactac acaatcaaag acgccctcgg cttacttctc
    14881 ttcattctct ccttaatgac attaacacta ttctcaccag acctcctagg cgacccagac
    14941 aattataccc tagccaaccc cttaaacacc cctccccaca tcaagcccga atgatatttc
    15001 ctattcgcct acacaattct ccgatccgtc cctaacaagc taggaggcgt ccttgcccta
    15061 ttactatcca tcctcatcct agcaataatc cccatcctcc atatatccaa acaacaaagc
    15121 ataatatttc gcccactaag ccaatcactt tattgactcc tagccgcaga cctcctcatt
    15181 ctaacctgaa tcggaggaca accagtaagc taccctttta ccatcattgg acaagtagca
    15241 tccgtactat acttcacaac aatcctaatc ctaataccaa ctatctccct aattgaaaac
    15301 aaaatactca aatgggcctg tccttgtagt ataaactaat acaccagtct tgtaaaccga
    15361 agatgaaaac ctttttccaa ggacaaatca gagaaaaagt ctttaactcc accattagca
    15421 cccaaagcta agattctaat ttaaactatt ctctgttctt tcatggggaa gcagatttgg
    15481 gtaccaccca agtattgact cacccatcaa caaccgctat gtatttcgta cattactgcc
    15541 agccaccatg aatattgcac ggtaccataa atacttgacc acctgtagta cataaaaacc
    15601 caatccacat caaaaccccc tccccatgct tacaagcaag tacagcaatc aaccctcaac
    15661 tatcacacat caactgcaac tccaaagcca cccctcaccc actaggatac caacaaacct
    15721 acctatcctt aacagtacat agtacataaa gccatttacc gtacatagca cattacagtc
    15781 aaatcccttc tcgtccccat ggatgacccc cctcagatag gggtcccttg accaccatcc
    15841 tccgtgaaat caatatcccg cacaagagtg ctactctcct cgctccgggc ccataacact
    15901 tgggggtagc taaagtgaac tgtatccgac atctggttcc tacttcaggg ccataaagcc
    15961 taaatagccc acacgttccc cttaaataag acatcacgat ggatcacagg tctatcaccc
    16021 tattaaccac tcacgggagc tctccatgca tttggtattt tcgtctgggg ggtgtgcacg
    16081 cgatagcatt gcgagacgct ggagccggag caccctatgt cgcagtatct gtctttgatt
    16141 cctgcctcat cctattattt atcgcaccta cgttcaatat tacaggcgaa catacttact
    16201 aaagtgtgtt aattaattaa tgcttgtagg acataataat aacaattgaa tgtctgtaca
    16261 gccgctttcc acacagacat cataacaaaa aatttccacc aaaccccccc tccccccgct
    16321 tctggccaca gcacttaaac acatctctgc caaaccccaa aaacaaagaa ccctaacacc
    16381 agcctaacca gatttcaaat tttatctttt ggcggtatgc acttttaaca gtcacccccc
    16441 aactaacaca ttattttccc ctcccactcc catactacta atctcatcaa tacaaccccc
    16501 gcccatccta cccagcacac acacaccgct gctaacccca taccccgaac caaccaaacc
    16561 ccaaagacac ccccccccac a
//
</code></pre> 
</details>







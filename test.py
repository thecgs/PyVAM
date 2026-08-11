import pyvam
import  matplotlib.pyplot as plt

#exmaple1
print("example1")
fig, ax = pyvam.draw_circos_MT("NC_012920.1", colors="OGDRAW", output="./doc/Fig.1.png", dpi=72)

#example2
print("example2")
fig, ax = pyvam.draw_circos_MT("MK804157", colors="OGDRAW", output="./doc/Fig.2.png", dpi=72)

#example3
print("example4")
fig, axs = plt.subplots(1, 3, figsize=(20, 20/3), subplot_kw={'projection':'polar'})
plt.subplots_adjust(wspace=0.3)

pyvam.draw_circos_MT(file="OP311642", colors="MitoFish", radius=20, gene_label_size=5, gene_label_inner=False, show_info=True, show_legend=False, axes=axs[0], info_fontsize=6)
pyvam.draw_circos_MT(file="OP289102", colors="MitoFish", radius=20, gene_label_size=5, gene_label_inner=False, show_info=True, show_legend=False, axes=axs[1], info_fontsize=6)
pyvam.draw_circos_MT(file="OP311642", colors="MitoFish", radius=20, gene_label_size=5, gene_label_inner=False, show_info=False,show_legend=False, show_GC_circos=False, axes=axs[2])
pyvam.draw_circos_MT(file="OP289102", colors="MitoFish", radius=12, gene_label_size=4, show_gene_label=True,   show_info=False,show_legend=False, show_GC_circos=False, axes=axs[2])

pyvam.add_tag(axs=axs, by_row=True)

axs[2].text(0.5, 0.5, s="Outer circle: Meghimatium pictum\nInner circle: Succinea arundinetorum", size=6, ha='center', va='center', style='italic')

fig.savefig("./doc/Fig.3.png", bbox_inches='tight', dpi=72)

#example4
print("example4")
pyvam.draw_linear_MT_nonproportional(files=["MK804148", "MK804158","MK804149", "MK804157"], 
                      start='ND1', add_id=True, dpi=72, force_reoriented=True,
                      output="./doc/Fig.4.png")
#example5
print("example5")
pyvam.draw_linear_MT(files=["MK804148", "MK804158","MK804149", "MK804157"], 
                      start='ND1', add_id=True, dpi=72, force_reoriented=True,
                      output="./doc/Fig.5.png")
#example6
print("example6")
pyvam.draw_linear_MT_nonproportional_interactive(files=["MK804148", "MK804158", "MK804149", "MK804157"],
                                                  start='COX1', add_id=True, force_reoriented=True,output="./doc/Fig.6.html")

pyvam.draw_linear_MT_interactive(files=["MK804148", "MK804158","MK804149", "MK804157"], 
                                  start='COX1', add_id=True, force_reoriented=True,
                                  output="./doc/Fig.7.html")

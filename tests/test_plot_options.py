from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

import pyvam


def write_minimal_record(path):
    record = SeqRecord(Seq("A" * 30), id="TEST", name="TEST")
    record.annotations["molecule_type"] = "DNA"
    record.features = [
        SeqFeature(
            SimpleLocation(0, 30, strand=1),
            type="source",
            qualifiers={"organism": ["Test species"]},
        ),
        SeqFeature(
            SimpleLocation(0, 27, strand=1),
            type="CDS",
            qualifiers={
                "gene": ["ND1"],
                "product": ["NADH dehydrogenase subunit 1"],
            },
        ),
    ]
    SeqIO.write(record, path, "genbank")


def test_interactive_custom_colors_and_label_color(tmp_path):
    input_path = tmp_path / "record.gbk"
    write_minimal_record(input_path)

    figure = pyvam.draw_linear_MT_nonproportional_interactive(
        input_path,
        colors={"source": "gray", "Other genes": "gray"},
        gene_label_color="red",
        show_legend=False,
        default_topology="linear",
    )

    annotation_colors = {
        annotation.font.color
        for annotation in figure.layout.annotations
        if annotation.font and annotation.font.color
    }
    assert annotation_colors == {"red"}


def test_static_api_accepts_default_topology(tmp_path):
    input_path = tmp_path / "record.gbk"
    write_minimal_record(input_path)

    figure, _ = pyvam.draw_circos_MT(
        input_path,
        default_topology="linear",
        show_legend=False,
        show_GC_circos=False,
    )
    assert figure is not None

"""Chloroplast detection, nomenclature, palette selection and plotting."""
import json
import re
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from matplotlib.colors import to_rgba

import pyvam
from pyvam.config import (CommonNamesDict_cp,
                          MTColors, cp_gene_color)
from pyvam.parserGB import ChloroplastWarning, _CP_NAMES, search_name


def cp_record(tmp_path, organelle="plastid:chloroplast"):
    record = SeqRecord(Seq("ATG" * 400), id="CP_TEST", name="CP_TEST")
    record.annotations.update(molecule_type="DNA", topology="circular")
    qualifiers = {"organism": ["Test species"]}
    if organelle:
        qualifiers["organelle"] = [organelle]
    record.features = [SeqFeature(SimpleLocation(0, 1200), type="source", qualifiers=qualifiers)]
    for i, (name, kind) in enumerate([
        ("photosystem I P700 apoprotein A1", "CDS"),
        ("rbcL", "CDS"), ("ndhA", "CDS"), ("atpB", "CDS"),
        ("16S ribosomal RNA", "rRNA"), ("tRNA-His", "tRNA"),
        ("orf999", "CDS"), ("unclassified protein", "CDS"),
        ("rrn23S", "gene"), ("trnH-GUG", "gene"),
    ]):
        record.features.append(SeqFeature(SimpleLocation(10 + i * 100, 70 + i * 100, strand=1),
                                          type=kind, qualifiers={"product": [name]}))
    path = tmp_path / "chloroplast.gb"
    SeqIO.write(record, path, "genbank")
    return path


def test_cp_dictionary_count_and_standard_identifiers():
    assert len(CommonNamesDict_cp) == 2167
    assert len(set(CommonNamesDict_cp.values())) == 102
    assert all(k == k.upper() for k in CommonNamesDict_cp)
    for canonical in set(CommonNamesDict_cp.values()):
        assert search_name(canonical, _CP_NAMES).casefold() == canonical.casefold()
    assert search_name("ATP synthase alpha subunit", _CP_NAMES) == "atpA"
    assert search_name("ATP synthase alpha subunit") != "atpA"


@pytest.mark.parametrize("organelle", ["chloroplast", "plastid:chloroplast", "Chloroplast"])
def test_cp_names_colors_and_rna_types(tmp_path, organelle):
    path = cp_record(tmp_path, organelle)
    with pytest.warns(ChloroplastWarning, match="CommonNamesDict_cp.*OGDRAW_CP"):
        features = pyvam.get_features(path)
    by_name = {f.name: f for f in features}
    assert features[0].organelle == "chloroplast"
    assert features[0].color_theme == "OGDRAW_CP"
    for name, color in [("psaA", "#00662C"), ("rbcL", "#1FA12D"), ("ndhA", "#FFEC00"),
                        ("atpB", "#97BE0D"), ("rrn16S", "#E2001A"), ("tRNA-His", "#162983"),
                        ("ORF999", "#57B9A8"), ("UNCLASSIFIED PROTEIN", "#AB259D")]:
        assert by_name[name].color == color
    assert by_name["rrn16S"].type == "rRNA"
    assert by_name["rrn23S"].type == "rRNA"
    assert len([f for f in features if f.name == "tRNA-His" and f.type == "tRNA"]) == 2


@pytest.mark.parametrize("organelle", [None, "mitochondrion", "plastid:apicoplast"])
def test_cp_is_not_inferred_from_gene_names(tmp_path, organelle):
    path = cp_record(tmp_path, organelle)
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        features = pyvam.get_features(path)
    assert not [w for w in captured if issubclass(w.category, ChloroplastWarning)]
    assert features[0].color_theme == "MITOFISH"


@pytest.mark.parametrize("colors,expected", [
    ("OGDRAW", "#808080ff"), ("AUTO", "#1FA12D"),
    ({"rbcL": "red", "Other genes": "gray"}, "red"),
])
def test_cp_explicit_palette_is_preserved(tmp_path, colors, expected):
    path = cp_record(tmp_path)
    with pytest.warns(ChloroplastWarning):
        features = pyvam.get_features(path, colors=colors)
    rbc = next(f for f in features if f.name == "rbcL")
    if isinstance(colors, str) and colors == "OGDRAW":
        expected = MTColors[colors]["Other genes"]
    assert rbc.color == expected


@pytest.mark.parametrize("name,kind,rgb", [
    ("psaA", "CDS", (0, 102, 44)), ("psbA", "CDS", (50, 137, 37)),
    ("petA", "CDS", (121, 156, 19)), ("atpA", "CDS", (151, 190, 13)),
    ("ndhA", "CDS", (255, 236, 0)), ("rbcL", "CDS", (31, 161, 45)),
    ("pbf1", "CDS", (254, 254, 200)), ("pafI", "CDS", (254, 254, 200)),
    ("pafII", "CDS", (254, 254, 200)), ("rpoA", "CDS", (189, 18, 32)),
    ("rps2", "CDS", (219, 170, 115)), ("rpl2", "CDS", (158, 119, 66)),
    ("unknown", "tRNA", (22, 41, 131)), ("unknown", "rRNA", (226, 0, 26)),
    ("clpP1", "CDS", (233, 93, 15)), ("matK", "CDS", (233, 93, 15)),
    ("ycf1", "CDS", (218, 218, 218)), ("orf1", "CDS", (87, 185, 168)),
    ("unknown", "rep_origin", (255, 128, 128)),
    ("unknown", "operon", (255, 0, 255)), ("unknown", "prim_transcript", (255, 0, 255)),
    ("unknown", "intron", (255, 255, 255)), ("unknown", "CDS", (171, 37, 157)),
])
def test_cp_palette_matches_supplied_ogdraw_rgb(name, kind, rgb):
    assert cp_gene_color(name, kind) == "#" + "".join(f"{c:02X}" for c in rgb)


@pytest.mark.parametrize("draw", [pyvam.draw_circos_MT, pyvam.draw_linear_MT,
                                  pyvam.draw_linear_MT_nonproportional])
def test_cp_static_default_palette_and_legend(tmp_path, draw):
    path = cp_record(tmp_path)
    options = {"show_xaxis": False} if draw == pyvam.draw_linear_MT else {}
    with pytest.warns(ChloroplastWarning):
        figure, axes = draw(path, **options)
    try:
        legend = axes.get_legend()
        labels = [text.get_text() for text in legend.get_texts()]
        assert "photosystem I" in labels
        assert "ribosomal RNAs" in labels
        assert "Complex IV (Cytochrome c oxidase)" not in labels
        assert any(p.get_facecolor() == to_rgba("#00662C") for p in axes.patches)
    finally:
        plt.close(figure)


@pytest.mark.parametrize("draw", [pyvam.draw_linear_MT_interactive,
                                  pyvam.draw_linear_MT_nonproportional_interactive])
def test_cp_interactive_legend_groups(tmp_path, draw):
    path = cp_record(tmp_path)
    with pytest.warns(ChloroplastWarning):
        figure = draw(path)
    assert any(trace.legendgroup == "photosystem I" and trace.fillcolor == "#00662C"
               for trace in figure.data)
    assert any(trace.legendgroup == "ribosomal RNAs" for trace in figure.data)


def test_cp_gene_cds_and_rna_placeholders_are_deduplicated(tmp_path):
    path = cp_record(tmp_path)
    record = SeqIO.read(path, "genbank")
    record.features = [record.features[0]]
    for kind, product, start in [("CDS", "rbcL", 10), ("rRNA", "16S ribosomal RNA", 110)]:
        location = SimpleLocation(start, start + 60, strand=1)
        record.features.extend([
            SeqFeature(location, type="gene", qualifiers={"gene": [product]}),
            SeqFeature(location, type=kind, qualifiers={"product": [product]}),
        ])
    SeqIO.write(record, path, "genbank")
    with pytest.warns(ChloroplastWarning):
        features = pyvam.get_features(path)
    assert [(f.name, f.type) for f in features[1:]] == [("rbcL", "CDS"), ("rrn16S", "rRNA")]


@pytest.mark.parametrize("view", ["circular", "linear", "order"])
def test_web_runtime_selected_cp_palette_and_warning(tmp_path, view):
    path = cp_record(tmp_path)
    source = (Path(__file__).parents[1] / "docs" / "app.js").read_text()
    runtime = source.split("const runtime = String.raw`", 1)[1].split("`;", 1)[0]
    runtime = runtime.replace("/home/pyvam-output", str(tmp_path / "output"))
    ns = {}
    exec(runtime, ns)
    options = json.loads(ns["pyvam_options"]())
    assert "OGDRAW_CP" in options["themes"] and "AUTO" not in options["themes"]
    assert "rbcL" not in options["markers"]
    assert "rbcL" in options["cp_markers"]
    assert set(options["markers"]).issubset(options["cp_markers"])
    assert options["class_colours"]["OGDRAW_CP"]["photosystem I"] == "#00662C"
    web = {key: None for key in re.findall(r'web\["([^\"]+)"\]', runtime)}
    web.update(dict.fromkeys(["gene_label_color", "species_label_size", "species_label_color", "height",
                              "radius", "info_fontsize", "GC_circos_height", "GC_circos_bin", "GC_circos_step",
                              "xaxisfontsize", "hspace", "subplot_height_cex"]))
    web.update(isfilename2species=False, abbr=False, add_id=True, remove_ncr=False,
               show_legend=True, force_reoriented=False, direction=-1, show_gc_circos=False,
               show_info=False, gene_label_inner=False, tidyname=False, show_xaxis=False)
    args = [json.dumps([str(path)]), view, "OGDRAW_CP", "", True, "{}", json.dumps(web)]
    assert ns["pyvam_render"](*args)
    assert any("detected chloroplast" in message for message in ns["_web_warnings"])
    legend = ns["_web_figure"].axes[0].get_legend()
    assert "photosystem I" in [text.get_text() for text in legend.get_texts()]
    plt.close("all")


@pytest.mark.parametrize("draw", [pyvam.draw_linear_MT, pyvam.draw_linear_MT_nonproportional,
                                  pyvam.draw_linear_MT_interactive,
                                  pyvam.draw_linear_MT_nonproportional_interactive])
def test_mixed_cp_and_mito_defaults_have_correct_legends(tmp_path, draw):
    cp_path = cp_record(tmp_path)
    record = SeqRecord(Seq("A" * 1200), id="MITO", name="MITO")
    record.annotations.update(molecule_type="DNA", topology="circular")
    record.features = [
        SeqFeature(SimpleLocation(0, 1200), type="source",
                   qualifiers={"organism": ["Mito species"], "organelle": ["mitochondrion"]}),
        SeqFeature(SimpleLocation(10, 70, strand=1), type="CDS", qualifiers={"gene": ["ND1"]}),
        SeqFeature(SimpleLocation(100, 160, strand=1), type="CDS", qualifiers={"gene": ["ATPase6"]}),
    ]
    mito_path = tmp_path / "mito.gb"
    SeqIO.write(record, mito_path, "genbank")
    options = {"show_xaxis": False} if draw == pyvam.draw_linear_MT else {}
    with pytest.warns(ChloroplastWarning):
        result = draw([cp_path, mito_path], **options)
    if draw in (pyvam.draw_linear_MT_interactive, pyvam.draw_linear_MT_nonproportional_interactive):
        entries = [t for t in result.data if t.showlegend]
        assert any(t.legendgroup == "photosystem I" for t in entries)
        assert any(t.legendgroup == "Complex I (NADH dehydrogenase)" for t in entries)
        assert len({t.fillcolor for t in entries if t.legendgroup == "ATP synthase"}) == 2
        return
    figure, axes = result
    try:
        legend_axis = figure.axes[-1] if draw == pyvam.draw_linear_MT else axes
        labels = [text.get_text() for text in legend_axis.get_legend().get_texts()]
        assert "photosystem I" in labels
        assert "Complex I (NADH dehydrogenase)" in labels
    finally:
        plt.close(figure)


def test_tidy_cp_retains_chloroplast_identity(tmp_path):
    path = cp_record(tmp_path)
    output = tmp_path / "tidy.gb"
    with pytest.warns(ChloroplastWarning):
        pyvam.tidy_genbank(path, output=output, table=11)
    record = SeqIO.read(output, "genbank")
    assert record.features[0].qualifiers["organelle"] == ["chloroplast"]
    with pytest.warns(ChloroplastWarning):
        features = pyvam.get_features(output)
    assert any(f.name == "rbcL" and f.color == "#1FA12D" for f in features)


def test_shared_cp_theme_recognizes_mt_genes_on_cp_record(tmp_path):
    path = cp_record(tmp_path)
    record = SeqIO.read(path, "genbank")
    for i, name in enumerate(["ND1", "COX1", "ATPase6", "Cytb"]):
        record.features.append(SeqFeature(SimpleLocation(1010 + i * 30, 1030 + i * 30, strand=1),
                                          type="CDS", qualifiers={"gene": [name]}))
    SeqIO.write(record, path, "genbank")
    with pytest.warns(ChloroplastWarning):
        features = pyvam.get_features(path, colors="OGDRAW_CP")
    by_name = {f.name: f for f in features}
    for name in ["ND1", "COX1", "ATPase6", "Cytb"]:
        assert name in by_name
        assert to_rgba(by_name[name].color) == to_rgba(MTColors["OGDRAW"][name])
    assert by_name["rbcL"].color == "#1FA12D"


@pytest.mark.parametrize("draw", [pyvam.draw_circos_MT, pyvam.draw_linear_MT,
                                  pyvam.draw_linear_MT_nonproportional])
def test_cp_legends_hide_absent_categories_and_balance_columns(tmp_path, draw):
    path = cp_record(tmp_path)
    options = {"show_xaxis": False} if draw == pyvam.draw_linear_MT else {}
    with pytest.warns(ChloroplastWarning):
        figure, axes = draw(path, colors="OGDRAW_CP", **options)
    try:
        legend = axes.get_legend()
        labels = [text.get_text() for text in legend.get_texts()]
        assert "photosystem I" in labels
        assert "RubisCO large subunit" in labels
        assert "photosystem II" not in labels
        assert "introns" not in labels
        assert "RNA polymerase" not in labels
        assert "cytochrome c oxidase" not in labels
        assert "Reverse" not in labels
        assert 1 < legend._ncols < 12
        assert len(set(labels)) == len(labels)
    finally:
        plt.close(figure)


@pytest.mark.parametrize("draw", [pyvam.draw_circos_MT, pyvam.draw_linear_MT,
                                  pyvam.draw_linear_MT_nonproportional])
def test_mt_legends_hide_absent_categories_even_with_short_labels(tmp_path, draw):
    path = cp_record(tmp_path, "mitochondrion")
    record = SeqIO.read(path, "genbank")
    record.features = [record.features[0],
                       SeqFeature(SimpleLocation(10, 70, strand=1), type="CDS",
                                  qualifiers={"gene": ["ND1"]})]
    SeqIO.write(record, path, "genbank")
    options = {"show_xaxis": False, "tidyname": True} if draw == pyvam.draw_linear_MT else {}
    figure, axes = draw(path, colors="OGDRAW_CP", **options)
    try:
        labels = [text.get_text() for text in axes.get_legend().get_texts()]
        assert "NADH dehydrogenase" in labels
        assert "photosystem I" not in labels
        assert "cytochrome c oxidase" not in labels
        assert "Reverse" not in labels
        assert "other genes" not in labels
    finally:
        plt.close(figure)


@pytest.mark.parametrize("draw", [pyvam.draw_linear_MT_interactive,
                                  pyvam.draw_linear_MT_nonproportional_interactive])
def test_shared_interactive_theme_groups_only_present_genes(tmp_path, draw):
    path = cp_record(tmp_path, "mitochondrion")
    record = SeqIO.read(path, "genbank")
    record.features = [record.features[0],
                       SeqFeature(SimpleLocation(10, 70, strand=1), type="CDS",
                                  qualifiers={"gene": ["COX1"]})]
    SeqIO.write(record, path, "genbank")
    figure = draw(path, colors="OGDRAW_CP")
    shown = [trace for trace in figure.data if trace.showlegend]
    assert any(trace.legendgroup == "cytochrome c oxidase" for trace in shown)
    assert not any(trace.name == "Reverse" or trace.legendgroup == "photosystem I" for trace in shown)
    assert figure.layout.legend.y == 0
    assert figure.layout.margin.b == 50



def test_shared_theme_accepts_cp_genes_without_changing_mt_names(tmp_path):
    path = cp_record(tmp_path, "mitochondrion")
    record = SeqIO.read(path, "genbank")
    record.features.append(SeqFeature(SimpleLocation(1100, 1160, strand=1), type="rRNA",
                                      qualifiers={"product": ["16S ribosomal RNA"]}))
    SeqIO.write(record, path, "genbank")
    with warnings.catch_warnings(record=True) as captured:
        features = pyvam.get_features(path, colors="OGDRAW_CP")
    assert not [w for w in captured if issubclass(w.category, ChloroplastWarning)]
    assert any(f.name == "rbcL" and f.color == "#1FA12D" for f in features)
    assert any(f.name == "16S rRNA" and f.type == "rRNA" for f in features)


@pytest.mark.parametrize("draw,columns,location,anchor", [
    (pyvam.draw_circos_MT, 1, 10, (1, -0.1)),
    (pyvam.draw_linear_MT, 12, 1, (1, 0.5)),
    (pyvam.draw_linear_MT_nonproportional, 12, 1, (1, 0)),
])
def test_mt_keeps_original_default_legend_layout(tmp_path, draw, columns, location, anchor):
    path = cp_record(tmp_path, "mitochondrion")
    options = {"show_xaxis": False} if draw == pyvam.draw_linear_MT else {}
    figure, axes = draw(path, colors="OGDRAW_CP", **options)
    try:
        legend = axes.get_legend()
        assert legend._ncols == columns
        assert legend._loc == location
        expected = axes.transAxes.transform(anchor)
        actual = legend.get_bbox_to_anchor()
        assert (actual.x0, actual.y0) == pytest.approx(tuple(expected))
    finally:
        plt.close(figure)


@pytest.mark.parametrize("length", [20000, 154478, 400000, 1200000])
@pytest.mark.parametrize("interactive", [False, True])
def test_cp_linear_ticks_scale_with_genome_length(tmp_path, length, interactive):
    path = cp_record(tmp_path)
    record = SeqIO.read(path, "genbank")
    record.seq = Seq("A" * length)
    record.features[0].location = SimpleLocation(0, length)
    SeqIO.write(record, path, "genbank")
    with pytest.warns(ChloroplastWarning):
        result = (pyvam.draw_linear_MT_interactive(path) if interactive
                  else pyvam.draw_linear_MT(path, show_xaxis=True))
    if interactive:
        step = result.layout.xaxis.dtick
    else:
        figure, axes = result
        try:
            ticks = axes[-1].get_xticks()
            step = ticks[1] - ticks[0]
            assert ticks[0] == 0 and ticks[-1] >= length
            assert len(ticks) <= 12
        finally:
            plt.close(figure)
    assert 5 <= length / step <= 10
    if length == 154478:
        assert step == 20000


def test_mixed_linear_uses_longest_genome_for_shared_cp_ticks(tmp_path):
    cp_path = cp_record(tmp_path)
    cp = SeqIO.read(cp_path, "genbank")
    cp.seq = Seq("A" * 154478)
    cp.features[0].location = SimpleLocation(0, 154478)
    SeqIO.write(cp, cp_path, "genbank")
    mito = SeqRecord(Seq("A" * 16569), id="MT", name="MT")
    mito.annotations.update(molecule_type="DNA", topology="circular")
    mito.features = [SeqFeature(SimpleLocation(0, 16569), type="source",
                               qualifiers={"organism": ["MT species"], "organelle": ["mitochondrion"]}),
                     SeqFeature(SimpleLocation(10, 100, strand=1), type="CDS", qualifiers={"gene": ["ND1"]})]
    mt_path = tmp_path / "mt.gb"
    SeqIO.write(mito, mt_path, "genbank")
    with pytest.warns(ChloroplastWarning):
        figure, axes = pyvam.draw_linear_MT([mt_path, cp_path])
    try:
        ticks = axes[-1].get_xticks()
        assert ticks[1] - ticks[0] == 20000
        assert list(axes[0].get_xticks()) == list(ticks)
    finally:
        plt.close(figure)


def test_mt_linear_ticks_keep_original_behavior(tmp_path):
    path = cp_record(tmp_path, "mitochondrion")
    figure, axes = pyvam.draw_linear_MT(path)
    try:
        ticks = axes[-1].get_xticks()
        assert ticks[1] - ticks[0] == 1000
    finally:
        plt.close(figure)
    figure = pyvam.draw_linear_MT_interactive(path)
    assert figure.layout.xaxis.dtick is None

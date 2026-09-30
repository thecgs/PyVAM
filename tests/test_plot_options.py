from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
import math

import pyvam
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
import pytest
from pyvam.drawMT import get_GC_bar_param, stat_features


@pytest.mark.parametrize("interactive", [False, True])
@pytest.mark.parametrize("topology", ["linear", "circular"])
@pytest.mark.parametrize("force", [None, False, True])
def test_proportional_rotation_respects_topology(tmp_path, interactive, topology, force):
    input_path = tmp_path / "rotation.gbk"
    record = SeqRecord(Seq("A" * 60), id="TEST", name="TEST")
    record.annotations.update(molecule_type="DNA", topology=topology)
    record.features = [
        SeqFeature(SimpleLocation(0, 60, strand=1), type="source",
                   qualifiers={"organism": ["Test species"]}),
        SeqFeature(SimpleLocation(0, 9, strand=1), type="CDS", qualifiers={"gene": ["ND2"]}),
        SeqFeature(SimpleLocation(12, 21, strand=1), type="CDS", qualifiers={"gene": ["ND1"]}),
    ]
    SeqIO.write(record, input_path, "genbank")
    options = {} if force is None else {"force_reoriented": force}
    rotated = topology == "circular" or force is True
    expected = {"ND1": 4.5, "ND2": 52.5 + (200 if topology == "linear" else 0)} if rotated else {
        "ND2": 4.5, "ND1": 16.5,
    }
    if interactive:
        figure = pyvam.draw_linear_MT_interactive(input_path, start="ND1", **options)
        positions = {a.text: a.x for a in figure.layout.annotations if a.text in expected}
        assert positions == expected
    else:
        figure, axes = pyvam.draw_linear_MT(input_path, start="ND1", show_xaxis=False, **options)
        try:
            positions = {t.get_text(): t.get_position()[0] for t in axes.texts if t.get_text() in expected}
            assert positions == expected
        finally:
            plt.close(figure)


@pytest.mark.parametrize("interactive", [False, True])
@pytest.mark.parametrize("reverse_order", [False, True])
def test_nonproportional_keeps_colocated_compound_annotations(tmp_path, interactive, reverse_order):
    input_path = tmp_path / "joined.gbk"
    record = SeqRecord(Seq("A" * 30), id="TEST", name="TEST")
    record.annotations.update(molecule_type="DNA", topology="circular")
    location = CompoundLocation([SimpleLocation(3, 9, strand=1), SimpleLocation(15, 24, strand=1)])
    annotations = [
        SeqFeature(location, type="CDS", qualifiers={"gene": ["ORF1"]}),
        SeqFeature(location, type="tRNA", qualifiers={"product": ["tRNA-Phe"]}),
    ]
    if reverse_order:
        annotations.reverse()
    record.features = [SeqFeature(SimpleLocation(0, 30, strand=1), type="source",
                                  qualifiers={"organism": ["Test species"]})] + annotations
    SeqIO.write(record, input_path, "genbank")
    if interactive:
        figure = pyvam.draw_linear_MT_nonproportional_interactive(input_path)
        labels = [a.text for a in figure.layout.annotations]
        assert labels.count("ORF1") == labels.count("F") == 1
    else:
        figure, axes = pyvam.draw_linear_MT_nonproportional(input_path)
        try:
            labels = [t.get_text() for t in axes.texts]
            assert labels.count("ORF1") == labels.count("F") == 1
        finally:
            plt.close(figure)


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
    plt.close(figure)


@pytest.mark.parametrize("show_xaxis", [True, False])
@pytest.mark.parametrize("custom_colors", [True, False])
def test_linear_plot_legend_with_single_genome(tmp_path, show_xaxis, custom_colors):
    input_path = tmp_path / "record.gbk"
    write_minimal_record(input_path)
    colors = {"source": "black", "ND1": "red"} if custom_colors else "mitofish"

    figure, axes = pyvam.draw_linear_MT(
        input_path, colors=colors, show_xaxis=show_xaxis, default_topology="circular",
    )
    try:
        legend_ax = axes[-1] if show_xaxis else axes
        legend = legend_ax.get_legend()
        assert legend is not None
        labels = [label.get_text() for label in legend.get_texts()]
        assert labels
        if custom_colors:
            assert labels == ["ND1"]
            assert legend.get_patches()[0].get_facecolor() == to_rgba("red")
        figure.canvas.draw()
    finally:
        plt.close(figure)


@pytest.mark.filterwarnings("error:Attempting to set identical low and high xlims:UserWarning")
@pytest.mark.parametrize("topology", ["circular", "linear"])
@pytest.mark.parametrize("rows", [
    [(1,)],
    [(-1,)],
    [(1, 1)],
    [(-1, -1)],
    [(1, -1, 1)],
    [(-1, 1, -1)],
    [(1, 1, 1), (-1,)],
])
def test_nonproportional_limits_include_all_features(tmp_path, topology, rows):
    paths = []
    for index, strands in enumerate(rows):
        length = 12 * len(strands) + 6
        record = SeqRecord(Seq("A" * length), id=f"TEST{index}", name=f"TEST{index}")
        record.annotations.update(molecule_type="DNA", topology=topology)
        record.features = [SeqFeature(
            SimpleLocation(0, length, strand=1), type="source",
            qualifiers={"organism": ["Test species"]},
        )]
        for gene_index, strand in enumerate(strands):
            record.features.append(SeqFeature(
                SimpleLocation(12 * gene_index, 12 * gene_index + 9, strand=strand),
                type="CDS", qualifiers={"gene": [["ND1", "COX1", "ATP6"][gene_index]]},
            ))
        path = tmp_path / f"record{index}.gbk"
        SeqIO.write(record, path, "genbank")
        paths.append(path)

    figure, axes = pyvam.draw_linear_MT_nonproportional(paths, show_legend=False)
    try:
        assert len(axes.patches) == sum(len(strands) for strands in rows)
        assert len(axes.collections) == (len(rows) if topology == "linear" else 0)
        drawn_x = [float(x) for patch in axes.patches for x in patch.get_xy()[:, 0]]
        drawn_x.extend(float(x) for collection in axes.collections
                       for segment in collection.get_segments() for x in segment[:, 0])
        left, right = axes.get_xlim()
        assert left < min(drawn_x) < max(drawn_x) < right
        # The plot should fit the content, not rely on Matplotlib's expansion
        # of a zero-width axis to an arbitrary default range.
        assert right - left <= 1.2 * (max(drawn_x) - min(drawn_x))
        figure.canvas.draw()
    finally:
        plt.close(figure)


@pytest.mark.parametrize("draw", [pyvam.draw_circos_MT, pyvam.draw_linear_MT_nonproportional])
@pytest.mark.parametrize("colors", [
    {"ND1": "red"},
    {"Other genes": "gray", "ND1": "red"},
    {"ND1": "red", "Other genes": "gray"},
])
@pytest.mark.parametrize("has_other_gene", [False, True])
def test_custom_legend_preserves_named_genes(tmp_path, draw, colors, has_other_gene):
    input_path = tmp_path / "record.gbk"
    write_minimal_record(input_path)
    if has_other_gene:
        record = SeqIO.read(input_path, "genbank")
        record.features.append(SeqFeature(
            SimpleLocation(27, 30, strand=1), type="CDS",
            qualifiers={"gene": ["ORF1"]},
        ))
        SeqIO.write(record, input_path, "genbank")

    figure, axes = draw(input_path, colors=colors, default_topology="circular")
    try:
        legend = axes.get_legend()
        assert legend is not None
        patches = {patch.get_label(): patch for patch in legend.get_patches()}
        labels = [text.get_text() for text in legend.get_texts()]
        assert "ND1" in labels
        assert patches["ND1"].get_facecolor() == to_rgba("red")
        assert ("Other genes" in labels) == (has_other_gene and "Other genes" in colors)
    finally:
        plt.close(figure)


@pytest.mark.parametrize("draw", [pyvam.draw_circos_MT, pyvam.draw_linear_MT_nonproportional])
def test_plot_saves_figure_for_supplied_axes(tmp_path, draw):
    input_path, output_path = tmp_path / "record.gbk", tmp_path / "figure.png"
    write_minimal_record(input_path)
    kwargs = {"projection": "polar"} if draw == pyvam.draw_circos_MT else {}
    figure, axes = plt.subplots(subplot_kw=kwargs)
    try:
        draw(input_path, axes=axes, output=output_path, default_topology="circular")
        assert axes.patches
        assert output_path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    finally:
        plt.close(figure)


@pytest.mark.parametrize("by_row", [True, False])
def test_add_tag_handles_two_dimensional_axes(by_row):
    figure, axes = plt.subplots(2, 3)
    try:
        pyvam.add_tag(axes, by_row=by_row, labels=["A", "B", "C", "D", "E", "F"])
        actual = [[axis.texts[-1].get_text() for axis in row] for row in axes]
        assert actual == ([["A", "B", "C"], ["D", "E", "F"]] if by_row
                          else [["A", "C", "E"], ["B", "D", "F"]])
    finally:
        plt.close(figure)


def test_add_tag_validates_axes_and_labels():
    with pytest.raises(ValueError, match="axs must"):
        pyvam.add_tag()
    figure, axis = plt.subplots()
    try:
        pyvam.add_tag(axis, labels=["Z"])
        assert axis.texts[-1].get_text() == "Z"
        with pytest.raises(ValueError, match="labels must contain"):
            pyvam.add_tag([axis, axis], labels=["A"])
    finally:
        plt.close(figure)


@pytest.mark.parametrize("bin_size,step", [(0, 50), (-1, 50), (50, 0), (50, -1), (True, 50)])
def test_gc_parameters_require_positive_integers(tmp_path, bin_size, step):
    input_path = tmp_path / "gc.gbk"
    write_minimal_record(input_path)
    with pytest.raises(ValueError, match="GC (bin|step) must be a positive integer"):
        get_GC_bar_param(pyvam.get_features(input_path), bin=bin_size, step=step)


@pytest.mark.parametrize("colors", [[], 1, object()])
def test_invalid_colors_raise_clear_error(tmp_path, colors):
    input_path = tmp_path / "colors.gbk"
    write_minimal_record(input_path)
    with pytest.raises(TypeError, match="colors must be"):
        pyvam.get_features(input_path, colors=colors)


@pytest.mark.parametrize("interactive", [False, True])
@pytest.mark.parametrize("filtered", [False, True])
@pytest.mark.parametrize("with_other_record", [False, True])
def test_nonproportional_plot_keeps_empty_genome_row(tmp_path, interactive, filtered, with_other_record):
    input_path = tmp_path / "empty.gbk"
    write_minimal_record(input_path)
    record = SeqIO.read(input_path, "genbank")
    record.annotations["topology"] = "circular"
    record.features = record.features[:1]
    if filtered:
        record.features.append(SeqFeature(
            SimpleLocation(3, 12, strand=1), type="D-loop",
            qualifiers={"note": ["Control Region"]},
        ))
    SeqIO.write(record, input_path, "genbank")
    files = [input_path]
    if with_other_record:
        other_path = tmp_path / "nonempty.gbk"
        write_minimal_record(other_path)
        files.append(other_path)

    if interactive:
        figure = pyvam.draw_linear_MT_nonproportional_interactive(files, remove_NCR=filtered)
        assert len(figure.layout.yaxis.ticktext) == len(files)
        assert len([trace for trace in figure.data if trace.fill == "toself"]) == int(with_other_record)
    else:
        figure, axes = pyvam.draw_linear_MT_nonproportional(files, remove_NCR=filtered)
        try:
            assert len(axes.get_yticklabels()) == len(files)
            assert len(axes.patches) == int(with_other_record)
            assert axes.get_xlim()[0] < axes.get_xlim()[1]
            figure.canvas.draw()
        finally:
            plt.close(figure)


@pytest.mark.parametrize("draw", [pyvam.draw_linear_MT, pyvam.draw_linear_MT_nonproportional])
@pytest.mark.parametrize("reverse_order", [False, True])
def test_batch_legends_include_all_genomes(tmp_path, draw, reverse_order):
    standard, other = tmp_path / "standard.gbk", tmp_path / "other.gbk"
    write_minimal_record(standard)
    write_minimal_record(other)
    record = SeqIO.read(other, "genbank")
    record.annotations["topology"] = "linear"
    record.features[1].qualifiers = {"gene": ["ORF1"]}
    SeqIO.write(record, other, "genbank")
    files = [other, standard] if reverse_order else [standard, other]

    figure, axes = draw(files, default_topology="circular")
    try:
        legend_ax = axes[-1] if draw == pyvam.draw_linear_MT else axes
        labels = [text.get_text() for text in legend_ax.get_legend().get_texts()]
        assert "Other genes" in labels
        assert "// Break" in labels
    finally:
        plt.close(figure)


@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("show_xaxis", [False, True])
@pytest.mark.parametrize("rows", [1, 2])
def test_static_proportional_gene_label_color(tmp_path, strand, show_xaxis, rows):
    input_path = tmp_path / "input.gbk"
    write_minimal_record(input_path)
    record = SeqIO.read(input_path, "genbank")
    record.features[1].location = SimpleLocation(0, 27, strand=strand)
    record.annotations["topology"] = "circular"
    SeqIO.write(record, input_path, "genbank")
    figure, axes = pyvam.draw_linear_MT(
        [input_path] * rows, gene_label_color="#cc2233", species_label_color="blue",
        show_xaxis=show_xaxis,
    )
    try:
        labels = [text for ax in figure.axes for text in ax.texts if text.get_text() == "ND1"]
        assert len(labels) == rows
        assert all(to_rgba(text.get_color()) == to_rgba("#cc2233") for text in labels)
        species = [text for ax in figure.axes for text in ax.texts if text.get_text() == "Test species"]
        assert len(species) == rows
        assert all(to_rgba(text.get_color()) == to_rgba("blue") for text in species)
    finally:
        plt.close(figure)


@pytest.mark.parametrize("topology", ["linear", "circular"])
def test_circos_gc_track_aligns_with_gene_track(tmp_path, topology):
    input_path = tmp_path / "input.gbk"
    record = SeqRecord(Seq("G" * 1000), id="TEST", name="TEST")
    record.annotations.update(molecule_type="DNA", topology=topology)
    record.features = [
        SeqFeature(SimpleLocation(0, 1000, strand=1), type="source",
                   qualifiers={"organism": ["Test species"]}),
        SeqFeature(SimpleLocation(500, 550, strand=1), type="CDS", qualifiers={"gene": ["ND1"]}),
    ]
    SeqIO.write(record, input_path, "genbank")
    figure, axes = pyvam.draw_circos_MT(input_path, GC_circos_bin=50, GC_circos_step=50)
    try:
        gene = axes.containers[0].patches[0]
        gc = axes.containers[-1].patches
        assert len(gc) == 20
        assert gc[10].get_x() == pytest.approx(gene.get_x())
        assert gc[10].get_width() == pytest.approx(gene.get_width())
        angle = 1.9 * math.pi if topology == "linear" else 2 * math.pi
        assert gc[-1].get_x() + gc[-1].get_width() == pytest.approx(angle)
        assert angle <= math.radians(axes.get_thetamax()) + 1e-12
        centers, values, widths = get_GC_bar_param(pyvam.get_features(input_path))
        assert centers[10] == pytest.approx(gene.get_x() + gene.get_width() / 2)
        assert widths[10] == pytest.approx(gene.get_width())
        assert values == [1.0] * 20
    finally:
        plt.close(figure)


@pytest.mark.parametrize("reverse_order", [False, True])
@pytest.mark.parametrize("same_name", [False, True])
@pytest.mark.parametrize("joined", [False, True])
def test_circos_counts_distinct_annotations(tmp_path, reverse_order, same_name, joined):
    input_path = tmp_path / "input.gbk"
    location = (CompoundLocation([SimpleLocation(3, 9, strand=1), SimpleLocation(15, 24, strand=1)])
                if joined else SimpleLocation(3, 24, strand=1))
    annotations = [SeqFeature(location, type="CDS", qualifiers={"gene": ["ORF1"], "locus_tag": ["LOC_A"]})]
    annotations.append(SeqFeature(location, type="CDS" if same_name else "tRNA", qualifiers={
        "gene": ["ORF1"] if same_name else ["tRNA-Phe"], "locus_tag": ["LOC_B"],
    }))
    if reverse_order:
        annotations.reverse()
    record = SeqRecord(Seq("A" * 30), id="TEST", name="TEST")
    record.annotations.update(molecule_type="DNA", topology="circular")
    record.features = [SeqFeature(SimpleLocation(0, 30, strand=1), type="source",
                                  qualifiers={"organism": ["Test species"]})] + annotations
    SeqIO.write(record, input_path, "genbank")
    features = pyvam.get_features(input_path)
    assert stat_features(features) == ((0, 0, 2) if same_name else (1, 0, 1))
    figure, axes = pyvam.draw_circos_MT(input_path, show_GC_circos=False)
    try:
        expected = "2 PCGs; 0 rRNAs; 0 tRNAs" if same_name else "1 PCGs; 0 rRNAs; 1 tRNAs"
        assert any(expected in text.get_text() for text in axes.texts)
    finally:
        plt.close(figure)
    # Both nonproportional consumers must retain separate loci as well.
    figure, axes = pyvam.draw_linear_MT_nonproportional(input_path)
    try:
        assert len(axes.patches) == 2
    finally:
        plt.close(figure)
    interactive = pyvam.draw_linear_MT_nonproportional_interactive(input_path)
    assert len(interactive.layout.annotations) == 2


def test_circos_returns_supplied_axes_and_empty_inputs_fail(tmp_path):
    input_path = tmp_path / "record.gbk"
    write_minimal_record(input_path)
    figure, axes = plt.subplots(subplot_kw={"projection": "polar"})
    try:
        assert pyvam.draw_circos_MT(input_path, axes=axes, show_GC_circos=False) is axes
    finally:
        plt.close(figure)

    with pytest.raises(ValueError, match="at least one input file"):
        pyvam.draw_linear_MT([])
    with pytest.raises(ValueError, match="at least one input file"):
        pyvam.draw_linear_MT_nonproportional([])
    with pytest.raises(ValueError, match="at least one input file"):
        pyvam.draw_linear_MT_interactive([])
    with pytest.raises(ValueError, match="at least one input file"):
        pyvam.draw_linear_MT_nonproportional_interactive([])

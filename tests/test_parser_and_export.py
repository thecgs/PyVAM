import builtins
import warnings
from io import StringIO

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import (
    BeforePosition, CompoundLocation, ExactPosition, SeqFeature,
    SimpleLocation, AfterPosition,
)
from Bio.SeqRecord import SeqRecord

from pyvam.parserGB import get_features, tidy_genbank
import pyvam.parserGB as parser_module


def write_record(path, sequence, features, topology=None):
    record = SeqRecord(Seq(sequence), id="TEST", name="TEST")
    record.annotations["molecule_type"] = "DNA"
    if topology is not None:
        record.annotations["topology"] = topology
    record.features = features
    SeqIO.write(record, path, "genbank")


def source_feature(length):
    return SeqFeature(
        SimpleLocation(0, length, strand=1),
        type="source",
        qualifiers={"organism": ["Test species"]},
    )


def cds_feature(location):
    return SeqFeature(
        location,
        type="CDS",
        qualifiers={
            "gene": ["ND1"],
            "product": ["NADH dehydrogenase subunit 1"],
        },
    )


def test_missing_topology_uses_configurable_fallback(tmp_path, caplog):
    input_path = tmp_path / "no-topology.gbk"
    write_record(
        input_path,
        "A" * 30,
        [source_feature(30), cds_feature(SimpleLocation(0, 27, strand=1))],
    )

    features = get_features(input_path, default_topology="LINEAR")

    assert features[0].topology == "linear"
    assert features[-1].name == "Gap"
    assert "does not declare a valid topology" in caplog.text


def test_invalid_default_topology_is_rejected(tmp_path):
    input_path = tmp_path / "record.gbk"
    write_record(input_path, "A" * 30, [source_feature(30)])

    with pytest.raises(ValueError, match="default_topology"):
        get_features(input_path, default_topology="diagonal")


def test_tidy_genbank_preserves_joined_cds_translation(tmp_path):
    input_path = tmp_path / "joined.gbk"
    output_path = tmp_path / "tidied.gbk"
    sequence = "ATGAAAAAA" + "C" * 6 + "ATGAAAAAA" + "C" * 6
    joined_location = CompoundLocation(
        [SimpleLocation(0, 9, strand=1), SimpleLocation(15, 24, strand=1)]
    )
    # A broad gene feature before the multipart CDS reproduces a common
    # GenBank layout and guards against the placeholder-suppression regression.
    features = [
        source_feature(len(sequence)),
        SeqFeature(SimpleLocation(0, 24, strand=1), type="gene", qualifiers={"gene": ["ND1"]}),
        cds_feature(joined_location),
    ]
    write_record(input_path, sequence, features, topology="circular")

    tidy_genbank(input_path, output=output_path, table=1)

    record = SeqIO.read(output_path, "genbank")
    cds = next(feature for feature in record.features if feature.type == "CDS")
    assert len(cds.location.parts) == 2
    assert cds.qualifiers["translation"] == ["MKKMKK"]


def test_tidy_genbank_honours_requested_translation_table(tmp_path):
    input_path = tmp_path / "table-one.gbk"
    output_path = tmp_path / "table-one-out.gbk"
    # GTG is a start codon in table 2 but not table 1; the output must retain
    # V when table 1 is selected.
    write_record(
        input_path,
        "GTGAAAAAATAA" + "AAA",
        [source_feature(15), cds_feature(SimpleLocation(0, 12, strand=1))],
        topology="circular",
    )

    tidy_genbank(input_path, output=output_path, table=1)

    record = SeqIO.read(output_path, "genbank")
    cds = next(feature for feature in record.features if feature.type == "CDS")
    assert cds.qualifiers["translation"] == ["VKK"]


@pytest.mark.parametrize("wraps_origin", [False, True])
def test_negative_join_roundtrip_preserves_sequence(tmp_path, wraps_origin):
    input_path, output_path = tmp_path / "input.gbk", tmp_path / "output.gbk"
    sequence = Seq("TTATTTCAT" + "C" * 6 + "GGGTTTCAT" + "C" * 6)
    parts = [SimpleLocation(15, 24, strand=-1), SimpleLocation(0, 9, strand=-1)]
    location = CompoundLocation(parts[::-1] if wraps_origin else parts)
    write_record(input_path, sequence, [source_feature(30), cds_feature(location)], "circular")

    tidy_genbank(input_path, output=output_path, table=1)

    record = SeqIO.read(output_path, "genbank")
    cds = next(f for f in record.features if f.type == "CDS")
    expected = location.extract(sequence)
    assert cds.location == location
    assert cds.extract(record.seq) == expected
    assert cds.qualifiers["translation"] == [str(expected.translate(table=1)).removesuffix("*")]


@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("origin", [0, 6, 12, 21])
def test_rotation_preserves_joined_cds_sequence(tmp_path, strand, origin):
    input_path, output_path = tmp_path / "input.gbk", tmp_path / "output.gbk"
    sequence = Seq("ATGAAAAAA" + "C" * 6 + "ATGAAAAAA" + "C" * 6)
    parts = [SimpleLocation(0, 9, strand=strand), SimpleLocation(15, 24, strand=strand)]
    location = CompoundLocation(parts if strand == 1 else parts[::-1])
    marker = SeqFeature(SimpleLocation(origin, origin + 2, strand=1),
                        type="tRNA", qualifiers={"product": ["tRNA-Phe"]})
    write_record(input_path, sequence,
                 [source_feature(30), cds_feature(location), marker], "circular")

    tidy_genbank(input_path, output=output_path, start="tRNA-Phe", table=1)

    record = SeqIO.read(output_path, "genbank")
    coding_features = [f for f in record.features if f.type == "CDS"]
    assert len(coding_features) == 1
    cds = coding_features[0]
    assert record.seq == sequence[origin:] + sequence[:origin]
    assert len(cds.location) == 18
    expected = location.extract(sequence)
    assert cds.extract(record.seq) == expected
    assert cds.qualifiers["translation"] == [str(expected.translate(table=1)).removesuffix("*")]
    assert all(0 <= part.start < part.end <= len(sequence) for part in cds.location.parts)


@pytest.mark.parametrize("strand", [1, -1])
def test_rotation_splits_simple_cds_crossing_new_origin(tmp_path, strand):
    input_path, output_path = tmp_path / "input.gbk", tmp_path / "output.gbk"
    sequence = Seq("CCC" * 2 + "ATGAAACCCGGGTTTTAA" + "CCC" * 2)
    location = SimpleLocation(6, 24, strand=strand)
    marker = SeqFeature(SimpleLocation(12, 14, strand=1), type="tRNA",
                        qualifiers={"product": ["tRNA-Phe"]})
    write_record(input_path, sequence,
                 [source_feature(30), cds_feature(location), marker], "circular")

    tidy_genbank(input_path, output=output_path, start="tRNA-Phe", table=1)

    record = SeqIO.read(output_path, "genbank")
    cds = next(f for f in record.features if f.type == "CDS")
    assert len(cds.location.parts) == 2
    assert cds.extract(record.seq) == location.extract(sequence)


@pytest.mark.parametrize("codon_start", [2, 3])
@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("joined", [False, True])
def test_export_preserves_codon_start_after_rotation(tmp_path, codon_start, strand, joined):
    input_path, output_path = tmp_path / "input.gbk", tmp_path / "output.gbk"
    coding = Seq("A" * (codon_start - 1) + "GTGAAATAA")
    genomic = coding if strand == 1 else coding.reverse_complement()
    if joined:
        sequence = Seq("CCC") + genomic[:5] + Seq("CCC") + genomic[5:] + Seq("CCC")
        parts = [SimpleLocation(3, 8, strand=strand),
                 SimpleLocation(11, 11 + len(genomic) - 5, strand=strand)]
        location = CompoundLocation(parts if strand == 1 else parts[::-1])
    else:
        sequence = Seq("CCC") + genomic + Seq("CCC")
        location = SimpleLocation(3, 3 + len(genomic), strand=strand)
    cds = cds_feature(location)
    cds.qualifiers["codon_start"] = [str(codon_start)]
    # The gene placeholder must not hide the following CDS reading frame.
    gene = SeqFeature(location, type="gene", qualifiers={"gene": ["ND1"]})
    write_record(input_path, sequence, [source_feature(len(sequence)), gene, cds], "circular")

    tidy_genbank(input_path, output=output_path, start="ND1", table=2)

    record = SeqIO.read(output_path, "genbank")
    result = next(f for f in record.features if f.type == "CDS")
    assert result.extract(record.seq) == coding
    assert result.qualifiers["codon_start"] == [str(codon_start)]
    # A partial CDS must not reinterpret the first complete GTG as a start codon.
    assert result.qualifiers["translation"] == ["VK"]


def test_export_to_stdout_retains_rna_and_control_region(tmp_path, capsys):
    input_path = tmp_path / "input.gbk"
    annotations = [source_feature(60)] + [
        SeqFeature(SimpleLocation(start, start + 9, strand=-1), type=kind,
                   qualifiers={"product": [name]})
        for start, kind, name in [(3, "tRNA", "tRNA-Phe"),
                                  (15, "rRNA", "12S rRNA"),
                                  (30, "D-loop", "Control Region")]
    ]
    write_record(input_path, "A" * 60, annotations, "linear")

    tidy_genbank(input_path)

    record = SeqIO.read(StringIO(capsys.readouterr().out), "genbank")
    assert record.seq == Seq("A" * 60)
    assert record.annotations["topology"] == "linear"
    assert record.annotations["organism"] == "Test species"
    assert [f.type for f in record.features] == ["source", "tRNA", "rRNA", "D-loop"]
    assert [f.location for f in record.features] == [f.location for f in annotations]


@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("layout", ["simple", "joined", "origin_join"])
@pytest.mark.parametrize("partial_end", [None, "five_prime", "three_prime"])
@pytest.mark.parametrize("rotate", [False, True])
def test_start_codon_conversion_requires_complete_five_prime_end(
    tmp_path, strand, layout, partial_end, rotate,
):
    input_path, output_path = tmp_path / "input.gbk", tmp_path / "output.gbk"
    coding = Seq("GTGAAATAA")
    if layout == "simple":
        spans = [(3, 12)]
    elif layout == "joined":
        spans = [(3, 6), (15, 21)]
    else:
        spans = [(21, 24), (0, 6)]
    if strand == -1:
        spans = spans[::-1]
    sequence = list("C" * 30)
    parts = []
    offset = 0
    for index, (start, end) in enumerate(spans):
        fragment = coding[offset:offset + end - start]
        sequence[start:end] = str(fragment if strand == 1 else fragment.reverse_complement())
        offset += end - start
        start, end = ExactPosition(start), ExactPosition(end)
        if partial_end == "five_prime" and index == 0:
            if strand == 1:
                start = BeforePosition(start)
            else:
                end = AfterPosition(end)
        if partial_end == "three_prime" and index == len(spans) - 1:
            if strand == 1:
                end = AfterPosition(end)
            else:
                start = BeforePosition(start)
        parts.append(SimpleLocation(start, end, strand=strand))
    location = parts[0] if len(parts) == 1 else CompoundLocation(parts)
    # Position 4 cuts a CDS part, exercising preservation of partial boundaries
    # when rotation splits a location around the new origin.
    marker = SeqFeature(SimpleLocation(4, 5, strand=1), type="tRNA",
                        qualifiers={"product": ["tRNA-Phe"]})
    write_record(input_path, "".join(sequence),
                 [source_feature(30), cds_feature(location), marker], "circular")

    tidy_genbank(input_path, output=output_path, table=2,
                 start="tRNA-Phe" if rotate else None)

    record = SeqIO.read(output_path, "genbank")
    cds = next(f for f in record.features if f.type == "CDS")
    assert cds.extract(record.seq) == coding
    assert cds.qualifiers["translation"] == ["VK" if partial_end == "five_prime" else "MK"]
    first_part = cds.location.parts[0]
    five_prime = first_part.start if strand == 1 else first_part.end
    assert (type(five_prime) is ExactPosition) == (partial_end != "five_prime")


@pytest.mark.parametrize("name", ["ND1", "ORF1"])
@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("joined", [False, True])
@pytest.mark.parametrize("gene_first", [False, True])
def test_cds_replaces_broad_gene_placeholder(tmp_path, name, strand, joined, gene_first):
    input_path, output_path = tmp_path / "input.gbk", tmp_path / "output.gbk"
    coding = Seq("AATGAAATAA")
    genomic = coding if strand == 1 else coding.reverse_complement()
    if joined:
        sequence = Seq("CCC") + genomic[:4] + Seq("CCCCCC") + genomic[4:] + Seq("CCC" * 4)
        parts = [SimpleLocation(3, 7, strand=strand), SimpleLocation(13, 19, strand=strand)]
        location = CompoundLocation(parts if strand == 1 else parts[::-1])
    else:
        sequence = Seq("CCC") + genomic + Seq("CCC" * 4)
        location = SimpleLocation(3, 13, strand=strand)
    # The broad gene uses an alias for ND1, exercising canonical matching.
    gene = SeqFeature(SimpleLocation(0, location.end + 3, strand=strand), type="gene",
                      qualifiers={"gene": ["NAD1" if name == "ND1" else name]})
    cds = SeqFeature(location, type="CDS", qualifiers={"gene": [name], "codon_start": ["2"]})
    annotations = [gene, cds] if gene_first else [cds, gene]
    write_record(input_path, sequence, [source_feature(len(sequence))] + annotations, "circular")

    features = get_features(input_path)
    parsed_cds = [feature for feature in features if feature.type == "CDS"]
    assert len(parsed_cds) == len(location.parts)
    assert all(f.name == name and f.codon_start == 2 for f in parsed_cds)
    assert all((f.join if f.join is not None else f.location) == location for f in parsed_cds)

    tidy_genbank(input_path, output=output_path, table=1)

    record = SeqIO.read(output_path, "genbank")
    result = [f for f in record.features if f.type == "CDS"]
    assert len(result) == 1
    assert result[0].location == location
    assert result[0].extract(record.seq) == coding
    assert result[0].qualifiers["codon_start"] == ["2"]
    assert result[0].qualifiers["translation"] == ["MK"]


@pytest.mark.parametrize("strand", [1, -1])
def test_exact_gene_does_not_hide_partial_cds_boundary(tmp_path, strand):
    input_path, output_path = tmp_path / "input.gbk", tmp_path / "output.gbk"
    coding = Seq("GTGAAATAA")
    sequence = Seq("CCC") + (coding if strand == 1 else coding.reverse_complement()) + Seq("CCC")
    location = SimpleLocation(BeforePosition(3) if strand == 1 else 3,
                              AfterPosition(12) if strand == -1 else 12, strand=strand)
    gene = SeqFeature(SimpleLocation(3, 12, strand=strand), type="gene",
                      qualifiers={"gene": ["ND1"]})
    write_record(input_path, sequence, [source_feature(15), gene, cds_feature(location)], "circular")

    tidy_genbank(input_path, output=output_path, table=2)

    record = SeqIO.read(output_path, "genbank")
    cds = next(f for f in record.features if f.type == "CDS")
    assert str(cds.location) == str(location)
    assert cds.qualifiers["translation"] == ["VK"]


@pytest.mark.parametrize("gene_location", [
    SimpleLocation(9, 15, strand=1),  # Within the CDS envelope, but in its gap.
    SimpleLocation(30, 39, strand=1),  # A separate copy of the same gene.
    SimpleLocation(3, 6, strand=-1),  # Overlaps a CDS part on the other strand.
])
@pytest.mark.parametrize("gene_first", [False, True])
def test_cds_keeps_unrelated_gene_only_annotations(tmp_path, gene_location, gene_first):
    input_path = tmp_path / "input.gbk"
    location = CompoundLocation([SimpleLocation(0, 9, strand=1),
                                 SimpleLocation(15, 24, strand=1)])
    gene = SeqFeature(gene_location, type="gene", qualifiers={"gene": ["ND1"]})
    cds = cds_feature(location)
    write_record(input_path, "A" * 45,
                 [source_feature(45)] + ([gene, cds] if gene_first else [cds, gene]), "circular")

    features = get_features(input_path)

    assert len(features) == 4
    assert any(f.location == gene_location and f.join is None for f in features[1:])
    assert sum(f.join == location for f in features[1:]) == 2


def test_gene_and_cds_match_by_locus_tag(tmp_path):
    input_path = tmp_path / "input.gbk"
    gene = SeqFeature(SimpleLocation(0, 24, strand=1), type="gene",
                      qualifiers={"gene": ["ORF1"], "locus_tag": ["TEST_001"]})
    cds = SeqFeature(SimpleLocation(3, 21, strand=1), type="CDS",
                     qualifiers={"product": ["hypothetical protein"], "locus_tag": ["TEST_001"]})
    write_record(input_path, "A" * 30, [source_feature(30), gene, cds], "circular")

    features = get_features(input_path)

    assert len(features) == 2
    assert features[1].location == cds.location


def test_gene_only_record_is_preserved(tmp_path):
    input_path = tmp_path / "input.gbk"
    gene = SeqFeature(SimpleLocation(3, 21, strand=1), type="gene", qualifiers={"gene": ["ORF1"]})
    write_record(input_path, "A" * 30, [source_feature(30), gene], "circular")

    features = get_features(input_path)

    assert len(features) == 2
    assert features[1].name == "ORF1"
    assert features[1].location == gene.location


@pytest.mark.parametrize("name", ["ND1", "ORF1"])
@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("joined", [False, True])
@pytest.mark.parametrize("gene_first", [False, True])
def test_unnamed_cds_inherits_matching_gene_name(tmp_path, name, strand, joined, gene_first):
    input_path, output_path = tmp_path / "input.gbk", tmp_path / "output.gbk"
    coding = Seq("AATGAAATAA")
    genomic = coding if strand == 1 else coding.reverse_complement()
    if joined:
        sequence = Seq("CCC") + genomic[:4] + Seq("CCCCCC") + genomic[4:] + Seq("CCC" * 4)
        parts = [SimpleLocation(3, 7, strand=strand), SimpleLocation(13, 19, strand=strand)]
        location = CompoundLocation(parts if strand == 1 else parts[::-1])
    else:
        sequence = Seq("CCC") + genomic + Seq("CCC" * 4)
        location = SimpleLocation(3, 13, strand=strand)
    gene = SeqFeature(SimpleLocation(0, location.end + 3, strand=strand), type="gene",
                      qualifiers={"gene": [name], "locus_tag": ["TEST_001"]})
    cds = SeqFeature(location, type="CDS", qualifiers={
        "locus_tag": ["TEST_001"], "protein_id": ["TEST_PROTEIN"],
        "codon_start": ["2"], "translation": ["MK"],
    })
    write_record(input_path, sequence,
                 [source_feature(len(sequence))] + ([gene, cds] if gene_first else [cds, gene]),
                 "circular")

    features = get_features(input_path)
    assert len(features) == 1 + len(location.parts)
    assert all(f.name == name and f.codon_start == 2 for f in features[1:])

    tidy_genbank(input_path, output=output_path, table=1)

    result = SeqIO.read(output_path, "genbank")
    coding_features = [f for f in result.features if f.type == "CDS"]
    assert len(coding_features) == 1
    exported = coding_features[0]
    assert exported.location == location
    assert exported.extract(result.seq) == coding
    assert exported.qualifiers["gene"] == [name]
    assert exported.qualifiers["codon_start"] == ["2"]
    assert exported.qualifiers["translation"] == ["MK"]


def test_unnamed_cds_does_not_take_name_from_unrelated_gene(tmp_path):
    input_path = tmp_path / "input.gbk"
    unrelated = SeqFeature(SimpleLocation(24, 33, strand=1), type="gene",
                           qualifiers={"gene": ["ND2"], "locus_tag": ["TEST_002"]})
    matching = SeqFeature(SimpleLocation(0, 21, strand=1), type="gene",
                          qualifiers={"gene": ["ND1"], "locus_tag": ["TEST_001"]})
    cds = SeqFeature(SimpleLocation(3, 21, strand=1), type="CDS",
                     qualifiers={"locus_tag": ["TEST_001"], "translation": ["KKKKKK"]})
    write_record(input_path, "A" * 45, [source_feature(45), unrelated, matching, cds], "circular")

    features = get_features(input_path)

    assert [(f.name, f.location) for f in features[1:]] == [
        ("ND1", cds.location), ("ND2", unrelated.location),
    ]


@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("annotation_types", [("CDS",), ("gene",), ("gene", "CDS")])
def test_full_length_gene_is_not_a_duplicate_of_source(tmp_path, strand, annotation_types):
    input_path, output_path = tmp_path / "input.gbk", tmp_path / "output.gbk"
    coding = Seq("ATG" + "AAA" * 8 + "TAA")
    location = SimpleLocation(0, len(coding), strand=strand)
    sequence = coding if strand == 1 else coding.reverse_complement()
    annotations = [SeqFeature(location, type=kind, qualifiers={"gene": ["ND1"]})
                   for kind in annotation_types]
    write_record(input_path, sequence, [source_feature(len(sequence))] + annotations, "circular")

    features = get_features(input_path)
    assert [(f.name, f.type) for f in features] == [("Test species", "source"), ("ND1", "CDS")]
    tidy_genbank(input_path, output=output_path, table=1)

    record = SeqIO.read(output_path, "genbank")
    cds = [f for f in record.features if f.type == "CDS"]
    assert len(cds) == 1
    assert cds[0].extract(record.seq) == coding
    assert cds[0].qualifiers["translation"] == ["M" + "K" * 8]


@pytest.mark.parametrize("reverse_order", [False, True])
def test_colocated_rna_and_cds_are_distinct_features(tmp_path, reverse_order):
    input_path = tmp_path / "input.gbk"
    location = SimpleLocation(0, 30, strand=1)
    annotations = [cds_feature(location), SeqFeature(location, type="tRNA",
                   qualifiers={"product": ["tRNA-Phe"]})]
    if reverse_order:
        annotations.reverse()
    # Repeating the CDS itself should still be deduplicated.
    annotations.append(cds_feature(location))
    write_record(input_path, "A" * 30, [source_feature(30)] + annotations, "circular")

    features = get_features(input_path)

    assert len(features) == 3
    assert {(f.name, f.type) for f in features[1:]} == {("ND1", "CDS"), ("tRNA-Phe", "tRNA")}


@pytest.mark.parametrize("name", ["ND1", "ORF1"])
@pytest.mark.parametrize("description", [
    {"product": ["hypothetical protein"]},
    {"note": ["predicted by annotation pipeline"]},
    {"product": ["hypothetical protein"], "note": ["predicted by annotation pipeline"]},
])
@pytest.mark.parametrize("gene_first", [False, True])
def test_cds_descriptions_do_not_replace_matching_gene_identifier(tmp_path, name, description, gene_first):
    input_path, output_path = tmp_path / "input.gbk", tmp_path / "output.gbk"
    gene = SeqFeature(SimpleLocation(3, 24, strand=1), type="gene",
                      qualifiers={"gene": [name], "locus_tag": ["TEST_001"]})
    location = SimpleLocation(6, 24, strand=1)
    cds = SeqFeature(location, type="CDS", qualifiers={
        "locus_tag": ["TEST_001"], "codon_start": ["1"], **description,
    })
    sequence = Seq("CCC" * 2 + "ATG" + "AAA" * 4 + "TAA" + "CCC" * 2)
    write_record(input_path, sequence,
                 [source_feature(30)] + ([gene, cds] if gene_first else [cds, gene]), "circular")

    features = get_features(input_path, start=name, colors={name: "red"})
    assert len(features) == 2
    assert features[1].name == name
    assert features[1].color == "red"
    assert features[1].location.start == 0
    assert features[0].mtgenome == sequence[6:] + sequence[:6]
    tidy_genbank(input_path, output=output_path, start=name, table=1)

    record = SeqIO.read(output_path, "genbank")
    result = next(f for f in record.features if f.type == "CDS")
    assert result.qualifiers["gene"] == [name]
    assert result.qualifiers["translation"] == ["MKKKK"]


def test_explicit_cds_gene_identifier_is_not_overwritten(tmp_path):
    input_path = tmp_path / "input.gbk"
    gene = SeqFeature(SimpleLocation(0, 24, strand=1), type="gene",
                      qualifiers={"gene": ["OLD_NAME"], "locus_tag": ["TEST_001"]})
    cds = SeqFeature(SimpleLocation(3, 21, strand=1), type="CDS", qualifiers={
        "gene": ["ORF1"], "locus_tag": ["TEST_001"], "product": ["hypothetical protein"],
        "note": ["predicted protein"],
    })
    write_record(input_path, "A" * 30, [source_feature(30), gene, cds], "circular")

    features = get_features(input_path)

    assert len(features) == 2
    assert features[1].name == "ORF1"


@pytest.mark.parametrize("kind", ["D-loop", "D_loop"])
@pytest.mark.parametrize("joined", [False, True])
def test_control_region_types_export_consistently(tmp_path, kind, joined):
    input_path, output_path = tmp_path / "input.gbk", tmp_path / "output.gbk"
    location = (CompoundLocation([SimpleLocation(21, 30, strand=1),
                                  SimpleLocation(0, 6, strand=1)]) if joined
                else SimpleLocation(3, 12, strand=1))
    control = SeqFeature(location, type=kind, qualifiers={"note": ["Control Region"]})
    write_record(input_path, "A" * 30, [source_feature(30), control], "circular")

    features = get_features(input_path)
    assert all(f.name == "D-loop" and f.type == "D-loop" for f in features[1:])
    tidy_genbank(input_path, output=output_path)

    result = SeqIO.read(output_path, "genbank")
    controls = [f for f in result.features if f.type == "D-loop"]
    assert len(controls) == 1
    assert controls[0].location == location
    assert controls[0].qualifiers["note"] == ["Control Region"]


@pytest.mark.parametrize("case", ["empty_file", "no_features", "no_source"])
def test_invalid_annotation_input_has_clear_error(tmp_path, case):
    input_path, output_path = tmp_path / "input.gbk", tmp_path / "output.gbk"
    if case == "empty_file":
        SeqIO.write([], input_path, "genbank")
    else:
        annotations = [] if case == "no_features" else [cds_feature(SimpleLocation(3, 21, strand=1))]
        write_record(input_path, "A" * 30, annotations, "circular")
    message = "source feature" if case == "no_source" else "feature annotations"

    with pytest.raises(ValueError, match=message):
        get_features(input_path)
    with pytest.raises(ValueError, match=message):
        tidy_genbank(input_path, output=output_path)
    assert not output_path.exists()


@pytest.mark.parametrize("kind,identifier,product,name", [
    ("tRNA", "MT-TF", "tRNA-Phe", "tRNA-Phe"),
    ("rRNA", "MT-RNR1", "12S ribosomal RNA", "12S rRNA"),
])
@pytest.mark.parametrize("joined", [False, True])
@pytest.mark.parametrize("known_gene", [False, True])
@pytest.mark.parametrize("gene_position", [None, "before", "after"])
def test_rna_product_fallback(tmp_path, kind, identifier, product, name, joined, known_gene, gene_position):
    input_path, output_path = tmp_path / "input.gbk", tmp_path / "output.gbk"
    location = (CompoundLocation([SimpleLocation(3, 9, strand=1),
                                  SimpleLocation(15, 24, strand=1)]) if joined
                else SimpleLocation(3, 24, strand=1))
    annotation = SeqFeature(location, type=kind, qualifiers={
        "gene": [name if known_gene else identifier],
        # A recognized gene name must still take precedence over the product.
        "product": ["tRNA-Val" if known_gene else product],
    })
    annotations = [annotation]
    if gene_position:
        gene = SeqFeature(location, type="gene", qualifiers={"gene": annotation.qualifiers["gene"]})
        annotations.insert(0 if gene_position == "before" else 1, gene)
    write_record(input_path, "A" * 30, [source_feature(30)] + annotations, "circular")

    features = get_features(input_path, colors={name: "red"})
    assert len(features) == (3 if joined else 2)
    assert all((f.name, f.type, f.color) == (name, kind, "red") for f in features[1:])
    tidy_genbank(input_path, output=output_path)
    exported = [f for f in SeqIO.read(output_path, "genbank").features if f.type == kind]
    assert len(exported) == 1
    assert exported[0].location == location
    assert exported[0].qualifiers["product"] == [product]


@pytest.mark.parametrize("reverse_order", [False, True])
@pytest.mark.parametrize("joined", [False, True])
@pytest.mark.parametrize("strand", [1, -1])
def test_colocated_orf_and_rna_keep_identity_and_translation(tmp_path, reverse_order, joined, strand):
    input_path, output_path = tmp_path / "input.gbk", tmp_path / "output.gbk"
    parts = [SimpleLocation(3, 7, strand=strand), SimpleLocation(13, 19, strand=strand)]
    if strand == -1:
        parts.reverse()
    location = CompoundLocation(parts) if joined else SimpleLocation(3, 13, strand=strand)
    coding = Seq("CATGAAATAA")  # codon_start=2 -> MK, with no intronic sequence.
    sequence = list("C" * 30)
    offset = 0
    for part in location.parts:
        fragment = coding[offset:offset + len(part)]
        if strand == -1:
            fragment = fragment.reverse_complement()
        sequence[int(part.start):int(part.end)] = str(fragment)
        offset += len(part)
    annotations = [
        SeqFeature(location, type="CDS", qualifiers={
            "gene": ["ORF1"], "locus_tag": ["LOC1"], "codon_start": ["2"],
        }),
        SeqFeature(location, type="tRNA", qualifiers={
            "product": ["tRNA-Phe"], "locus_tag": ["LOC2"],
        }),
    ]
    if reverse_order:
        annotations.reverse()
    write_record(input_path, "".join(sequence), [source_feature(30)] + annotations, "circular")

    features = get_features(input_path, colors={"ORF1": "red", "tRNA-Phe": "blue"})
    assert len(features) == (5 if joined else 3)
    assert {(f.name, f.type, f.color) for f in features[1:]} == {
        ("ORF1", "CDS", "red"), ("tRNA-Phe", "tRNA", "blue"),
    }
    assert all(f.codon_start == 2 for f in features if f.type == "CDS")
    tidy_genbank(input_path, output=output_path, table=1)
    record = SeqIO.read(output_path, "genbank")
    cds = [f for f in record.features if f.type == "CDS"]
    rna = [f for f in record.features if f.type == "tRNA"]
    assert len(cds) == len(rna) == 1
    assert cds[0].location == rna[0].location == location
    assert cds[0].qualifiers["gene"] == ["ORF1"]
    assert cds[0].qualifiers["translation"] == ["MK"]
    assert rna[0].qualifiers["product"] == ["tRNA-Phe"]


@pytest.mark.parametrize("reverse_order", [False, True])
def test_colocated_cds_keep_their_own_codon_start(tmp_path, reverse_order):
    input_path = tmp_path / "input.gbk"
    location = SimpleLocation(3, 24, strand=1)
    annotations = [
        SeqFeature(location, type="CDS", qualifiers={"gene": [name], "codon_start": [offset]})
        for name, offset in [("ORF1", "2"), ("ND1", "1")]
    ]
    if reverse_order:
        annotations.reverse()
    write_record(input_path, "A" * 30, [source_feature(30)] + annotations, "circular")

    features = get_features(input_path)

    assert len(features) == 3
    assert {(f.name, f.codon_start) for f in features[1:]} == {("ORF1", 2), ("ND1", 1)}


@pytest.mark.parametrize("remote", [False, True])
@pytest.mark.parametrize("case", [
    "success", "no_features", "invalid_codon_start", "parser_init_error", "parser_iteration_error",
])
def test_genbank_input_handle_is_closed(tmp_path, monkeypatch, remote, case):
    input_path = tmp_path / "input.gbk"
    cds = cds_feature(SimpleLocation(3, 24, strand=1))
    if case == "invalid_codon_start":
        cds.qualifiers["codon_start"] = ["4"]
    annotations = [] if case == "no_features" else [source_feature(30), cds]
    write_record(input_path, "A" * 30, annotations, "circular")
    opened = []

    def open_input(*args, **kwargs):
        handle = builtins.open(input_path)
        opened.append(handle)
        return handle

    if remote:
        # Simulate an NCBI response without making any network requests.
        monkeypatch.setattr(parser_module, "get_genbank_from_ncbi", open_input)
        file = "TEST_ACCESSION"
    else:
        monkeypatch.setattr(parser_module, "open", open_input, raising=False)
        original_parse = SeqIO.parse

        def parse_input(source, *args, **kwargs):
            # Track the old implicit-open path as well as the explicit handle.
            if source == input_path:
                source = open_input()
            return original_parse(source, *args, **kwargs)

        monkeypatch.setattr(parser_module.SeqIO, "parse", parse_input)
        file = input_path

    if case == "parser_init_error":
        def failing_parse(*args, **kwargs):
            raise ValueError("parser initialization failed")
        monkeypatch.setattr(parser_module.SeqIO, "parse", failing_parse)
    elif case == "parser_iteration_error":
        def failing_parse(*args, **kwargs):
            record = SeqRecord(Seq("A" * 30), id="TEST")
            record.annotations["topology"] = "circular"
            record.features = annotations
            yield record
            raise ValueError("parser iteration failed")
        monkeypatch.setattr(parser_module.SeqIO, "parse", failing_parse)

    try:
        if case == "success":
            assert len(get_features(file)) == 2
        else:
            messages = {
                "no_features": "no feature annotations",
                "invalid_codon_start": "CDS codon_start",
                "parser_init_error": "parser initialization failed",
                "parser_iteration_error": "parser iteration failed",
            }
            with pytest.raises(ValueError, match=messages[case]):
                get_features(file)
        assert len(opened) == 1
        assert opened[0].closed
    finally:
        # Also clean up when testing a regressed implementation.
        for handle in opened:
            handle.close()


@pytest.mark.parametrize("gene_first", [False, True])
@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("name", ["ND1", "ORF1"])
def test_nested_annotations_with_different_loci_remain_distinct(tmp_path, gene_first, strand, name):
    input_path, output_path = tmp_path / "input.gbk", tmp_path / "output.gbk"
    sequence = Seq("CCC" + "CATGAAATAA" + "CC" + "AAA" * 5)
    gene_location, cds_location = SimpleLocation(0, 15, strand=1), SimpleLocation(3, 13, strand=1)
    if strand == -1:
        sequence = sequence.reverse_complement()
        gene_location, cds_location = SimpleLocation(15, 30, strand=-1), SimpleLocation(17, 27, strand=-1)
    gene = SeqFeature(gene_location, type="gene", qualifiers={"gene": [name], "locus_tag": ["LOC_A"]})
    cds = SeqFeature(cds_location, type="CDS", qualifiers={
        "gene": [name], "locus_tag": ["LOC_B"], "codon_start": ["2"],
    })
    write_record(input_path, sequence,
                 [source_feature(30)] + ([gene, cds] if gene_first else [cds, gene]), "circular")

    features = get_features(input_path)
    assert len(features) == 3
    by_tag = {f.locus_tags[0]: f for f in features[1:]}
    assert by_tag["LOC_A"].location == gene_location
    assert by_tag["LOC_A"].original_type == "gene"
    assert by_tag["LOC_B"].location == cds_location
    assert by_tag["LOC_B"].original_type == "CDS"
    assert by_tag["LOC_B"].codon_start == 2
    tidy_genbank(input_path, output=output_path, table=1)
    record = SeqIO.read(output_path, "genbank")
    exported = {f.qualifiers["locus_tag"][0]: f for f in record.features if f.type == "CDS"}
    assert set(exported) == {"LOC_A", "LOC_B"}
    assert exported["LOC_A"].qualifiers["translation"] == ["PHEIT"]
    assert exported["LOC_B"].qualifiers["translation"] == ["MK"]
    assert exported["LOC_B"].location == cds_location


@pytest.mark.parametrize("reverse_order", [False, True])
@pytest.mark.parametrize("start", [None, "ND2"])
def test_same_named_compound_loci_survive_rotation_and_export(tmp_path, reverse_order, start):
    input_path, output_path = tmp_path / "input.gbk", tmp_path / "output.gbk"
    sequence = Seq("CCC" + "ATG" + "CCC" * 2 + "AAATAA" + "CCC" * 4)
    location = CompoundLocation([SimpleLocation(3, 6, strand=1), SimpleLocation(12, 18, strand=1)])
    annotations = [SeqFeature(location, type="CDS", qualifiers={
        "gene": ["ORF1"], "locus_tag": [tag],
    }) for tag in ("LOC_A", "LOC_B")]
    if reverse_order:
        annotations.reverse()
    # A repeated annotation must still be removed, unlike a distinct locus.
    annotations.append(annotations[0])
    annotations.append(SeqFeature(SimpleLocation(9, 12, strand=1), type="CDS", qualifiers={"gene": ["ND2"]}))
    write_record(input_path, sequence, [source_feature(30)] + annotations, "circular")

    features = get_features(input_path, start=start)
    orfs = [f for f in features if f.name == "ORF1"]
    assert len(orfs) == 4
    assert [f.locus_tags for f in orfs].count(("LOC_A",)) == 2
    assert [f.locus_tags for f in orfs].count(("LOC_B",)) == 2
    tidy_genbank(input_path, output=output_path, start=start, table=1)
    record = SeqIO.read(output_path, "genbank")
    exported = [f for f in record.features if f.type == "CDS" and f.qualifiers["gene"] == ["ORF1"]]
    assert len(exported) == 2
    assert {f.qualifiers["locus_tag"][0] for f in exported} == {"LOC_A", "LOC_B"}
    assert all(f.qualifiers["translation"] == ["MK"] for f in exported)
    assert all(f.extract(record.seq) == Seq("ATGAAATAA") for f in exported)
    assert len([f for f in get_features(output_path) if f.name == "ORF1"]) == 4


@pytest.mark.parametrize("reverse_order", [False, True])
def test_same_named_colocated_cds_keep_independent_frames(tmp_path, reverse_order):
    input_path = tmp_path / "input.gbk"
    location = SimpleLocation(3, 24, strand=1)
    annotations = [SeqFeature(location, type="CDS", qualifiers={
        "gene": ["ND1"], "locus_tag": [tag], "codon_start": [offset],
    }) for tag, offset in [("LOC_A", "1"), ("LOC_B", "2")]]
    if reverse_order:
        annotations.reverse()
    write_record(input_path, "A" * 30, [source_feature(30)] + annotations, "circular")
    features = get_features(input_path)
    assert len(features) == 3
    assert {(f.locus_tags, f.codon_start) for f in features[1:]} == {
        (("LOC_A",), 1), (("LOC_B",), 2),
    }


def test_source_annotation_order_and_partial_translation_are_safe(tmp_path):
    input_path, output_path = tmp_path / "input.gbk", tmp_path / "output.gbk"
    source = source_feature(5)
    cds = SeqFeature(SimpleLocation(0, 5, strand=1), type="CDS",
                     qualifiers={"gene": ["ND1"]})
    write_record(input_path, "ATGAA", [cds, source], "circular")

    assert get_features(input_path)[0].type == "source"
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("error")
        tidy_genbank(input_path, output=output_path, table=1)
    assert not caught
    exported = SeqIO.read(output_path, "genbank")
    assert next(f for f in exported.features if f.type == "CDS").qualifiers["translation"] == ["M"]


def test_alias_matching_does_not_rewrite_unrelated_prefixes(tmp_path):
    input_path = tmp_path / "input.gbk"
    write_record(input_path, "A" * 30, [source_feature(30),
                 SeqFeature(SimpleLocation(0, 9, strand=1), type="CDS",
                            qualifiers={"gene": ["COX10"]})], "circular")
    assert get_features(input_path)[1].name == "COX10"

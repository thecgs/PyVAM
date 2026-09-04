from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from pyvam.parserGB import get_features, tidy_genbank


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

    record = next(SeqIO.parse(output_path, "genbank"))
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

    record = next(SeqIO.parse(output_path, "genbank"))
    cds = next(feature for feature in record.features if feature.type == "CDS")
    assert cds.qualifiers["translation"] == ["VKK"]

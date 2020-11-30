"""
Generates a CLDF dataset for phoneme inventories from the "Journal of the IPA",
aggregated by Baird et al. (forth).
"""

import json
import unicodedata
from pathlib import Path
from unidecode import unidecode
import re

from pyglottolog import Glottolog
from pyclts import CLTS, models

from pycldf import Sources
from cldfbench import CLDFSpec
from cldfbench import Dataset as BaseDataset
from clldutils.misc import slug


def compute_id(text):
    """
    Returns a codepoint representation to an Unicode string.
    """

    unicode_repr = "".join(["u{0:0{1}X}".format(ord(char), 4) for char in text])

    label = slug(unidecode(text))

    return "%s_%s" % (label, unicode_repr)


def normalize_grapheme(text):
    """
    Apply simple, non-CLTS, normalization.
    """

    text = unicodedata.normalize("NFC", text)

    if text[0] == "(" and text[-1] == ")":
        text = text[1:-1]

    if text[0] == "[" and text[-1] == "]":
        text = text[1:-1]

    return text


def read_raw_source(filename):
    def _splitter(text):
        """
        Splits a list of phonemes as provided in the sourceself.

        We need to split by commas, provided they are not within parentheses (used to
        list allophones). This solution uses a negative-lookahead in regex.
        """
        items = [item.strip() for item in re.split(",\s*(?![^()]*\))", text)]
        return [item for item in items if item]

    # Holds the label to the current section of data
    section = None

    data = {
        "source": None,
        "language_name": None,
        "iso_code": None,
        "consonants": [],
        "vowels": [],
    }

    # Iterate over all lines
    with open(filename) as handler:
        for line in handler:
            # Clear line (including BOM) and skip empty data
            line = line.replace("\ufeff", "").strip()

            if not line:
                continue

            if line.startswith("#"):
                section = line[1:-1].strip()
            elif section == "Reference":
                data["source"] = line
            elif section == "Language":
                data["language_name"] = line
            elif section == "ISO Code":
                data["iso_code"] = line
            elif section == "Consonant Inventory":
                data["consonants"] += _splitter(line)
            elif section == "Vowel Inventory":
                data["vowels"] += _splitter(line)

    return data


class Dataset(BaseDataset):
    """
    CLDF dataset for inventories.
    """

    dir = Path(__file__).parent
    id = "jipa"

    def cldf_specs(self):  # A dataset must declare all CLDF sets it creates.
        return CLDFSpec(dir=self.cldf_dir, module="StructureDataset")

    def cmd_download(self, args):
        """
        Download files to the raw/ directory. You can use helpers methods of `self.raw_dir`, e.g.

        >>> self.raw_dir.download(url, fname)
        """
        pass

    def cmd_makecldf(self, args):
        """
        Convert the raw data to a CLDF dataset.

        >>> args.writer.objects['LanguageTable'].append(...)
        """

        # Instantiate Glottolog and CLTS
        # TODO: how to call CLTS?
        glottolog = Glottolog(args.glottolog.dir)
        clts_path = Path.home() / ".config" / "cldf" / "clts"
        clts_path = Path.home() / "src" / "INVENTORIES" / "clts"
        clts = CLTS(clts_path.absolute())
        clts_jipa = clts.transcriptiondata("jipa")

        # Add the bibliographic info
        sources = Sources.from_file(self.raw_dir / "sources.bib")
        args.writer.cldf.add_sources(*sources)

        # Add components
        args.writer.cldf.add_columns(
            "ValueTable",
            {"name": "Marginal", "datatype": "boolean"},
            "Catalog",
            "Contribution_ID",
        )

        args.writer.cldf.add_component("ParameterTable", "BIPA")
        args.writer.cldf.add_component(
            "LanguageTable", "Family_Glottocode", "Family_Name", "Glottolog_Name"
        )
        args.writer.cldf.add_table(
            "inventories.csv",
            "ID",
            "Name",
            "Contributor_ID",
            {
                "name": "Source",
                "propertyUrl": "http://cldf.clld.org/v1.0/terms.rdf#source",
                "separator": ";",
            },
            "URL",
            "Tones",
            primaryKey="ID",
        )

        # load language mapping and build inventory info
        languages = []
        for row in self.etc_dir.read_csv("languages.csv", dicts=True):
            if row["Glottocode"]:
                lang = glottolog.languoid(row["Glottocode"])
                update = {
                    "Family_Glottocode": lang.lineage[0][1] if lang.lineage else None,
                    "Family_Name": lang.lineage[0][0] if lang.lineage else None,
                    "Glottocode": row["Glottocode"],
                    "Latitude": lang.latitude,
                    "Longitude": lang.longitude,
                    "Macroarea": lang.macroareas[0].name if lang.macroareas else None,
                    "Glottolog_Name": lang.name,
                }
                row.update(update)

            languages.append(row)

        # Build source map from language
        source_map = {lang["ID"]: lang["Source"] for lang in languages}

        # Parse sources
        segments = []
        values = []
        counter = 1
        source_files = list(self.raw_dir.glob("*.txt"))
        for filename in source_files:
            contents = read_raw_source(filename)

            # TODO: Only keeping the main entry of allophones for now; we cannot just
            # strip parentheses and what is inside, as the notation (somewhat ambigous)
            # is also used for marginal sounds
            all_segments = []
            for segment in contents["consonants"] + contents["vowels"]:
                if segment[0] == "(" and segment[-1] == ")":
                    all_segments.append(segment)
                else:
                    stripped = re.sub(r"\([^)]*\)", "", segment)
                    all_segments.append(stripped)

            for segment in all_segments:
                # Obtain the corresponding BIPA grapheme, is possible
                normalized = normalize_grapheme(segment)

                # Due to the behavior of `.resolve_grapheme`, we need to attempt,
                # paying attention to raised exceptions, to convert in different ways
                sound = clts.bipa[
                    clts_jipa.grapheme_map.get(
                        segment, clts_jipa.grapheme_map.get(normalized, "")
                    )
                ]
                if isinstance(sound, models.UnknownSound):
                    sound = clts.bipa[normalized]

                if isinstance(sound, models.UnknownSound):
                    par_id = "UNK_" + compute_id(normalized)
                    bipa_grapheme = ""
                    desc = ""
                else:
                    par_id = "BIPA_" + compute_id(normalized)
                    bipa_grapheme = str(sound)
                    desc = sound.name
                segments.append((par_id, normalized, bipa_grapheme, desc))

                # TODO: do marginal
                lang_key = slug(contents["language_name"])
                values.append(
                    {
                        "ID": str(counter),
                        "Language_ID": lang_key,
                        "Marginal": False,
                        "Parameter_ID": par_id,
                        "Value": segment,
                        "Contribution_ID": lang_key,
                        "Source": [source_map[lang_key]],
                        "Catalog": "jipa",
                    }
                )
                counter += 1

        # Build segment data
        parameters = [
            {"ID": id, "Name": normalized, "BIPA": bipa_grapheme, "Description": desc}
            for id, normalized, bipa_grapheme, desc in set(segments)
        ]

        # Write data and validate
        inventories = []
        args.writer.write(
            **{
                "ValueTable": values,
                "LanguageTable": languages,
                "ParameterTable": parameters,
                "inventories.csv": inventories,
            }
        )

from urllib.request import urlopen
from Bio import Entrez as entrez
from Bio.Entrez.Parser import DictionaryElement, ListElement, StringElement
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from fnmatch import fnmatch
from io import BytesIO, TextIOWrapper
import openpyxl
from openpyxl import Workbook
from openpyxl.cell import Cell
from os import path
from pydantic import BaseModel, ValidationError
import re
from typing import NamedTuple, cast, Iterable, List, Literal, Optional, Sequence
from zipfile import ZipFile

from cbrextra.sequences.sequence import ConfigContext

from .data import SeqEntry
from .sequence import SeqEntryResult, SeqLoaderBase, SeqLoaderFactoryBase, SequenceLoadException

EXCEL_EXT = ".xlsx"
EXCEL_EXT_POS = -len(EXCEL_EXT)
ENTREZ_RETMAX = 1000
K_ID_LIST = "IdList"
K_ORGANISM_NAME = "Organism_Name"
K_ACCESSION = "Assembly_Accession"
K_DNA_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{" + K_ACCESSION + "}/download?include_annotation_type=GENOME_FASTA"

FIX_NAMES_RE = re.compile(r"\.(\w)")

class ExcelSheetSpec(BaseModel):
    name_pattern: str
    id_column: Optional[str]
    organism_column: str

class ExcelSpec(BaseModel):

    loader_type: Literal['excel']
    name_pattern: str
    sheets: List[ExcelSheetSpec]

class ExcelSpecConfig(NamedTuple):
    spec: ExcelSpec
    directory: str

class ExcelLoader(SeqLoaderBase):

    def __init__(self, configurations: Sequence[ExcelSpecConfig]):
        self.__configurations = configurations

    def __find_spec(self, file_path: str) -> Optional[ExcelSpec]:

        for configuration in self.__configurations:

            file_directory = path.dirname(file_path)
            file_name = path.basename(file_path)[0:EXCEL_EXT_POS]

            if file_directory == configuration.directory \
                and fnmatch(file_name, configuration.spec.name_pattern):
                return configuration.spec
            
        return None

    def load(self, file_path: str) -> Optional[Sequence[SeqEntryResult]]:

        if file_path[EXCEL_EXT_POS:].lower() != EXCEL_EXT:
            return None

        spec = self.__find_spec(file_path)

        if spec is None:
            return None

        return list(self.__load_from_file(file_path, spec))
    
    def __download_accession_gnomes(self, accession: str) -> Sequence[SeqRecord]:

        dna_zip_url = K_DNA_URL.format(**{K_ACCESSION: accession})
        with urlopen(dna_zip_url) as resp \
            , ZipFile(BytesIO(resp.read())) as archive:

            gnome_fasta_file = next(
                (file for file in archive.namelist() if file.endswith("genomic.fna")),
                None
            )

            if gnome_fasta_file is None:
                raise ValueError(f"The zip file at {dna_zip_url} has no DNA.")
            
            with archive.open(gnome_fasta_file, mode = "r") as fasta:

                return list(
                    SeqIO.parse(TextIOWrapper(fasta, encoding='utf-8'), format='fasta')
                )
    
    def __fix_name(self, name: str) -> str:

        # commas don't really mean much
        name = name.replace(",", " ")
        return FIX_NAMES_RE.sub(r". \g<1>", name)

    def __load_from_organism(self, organism: str, organism_id: Optional[str]) -> Iterable[SeqEntry]:

        organism = self.__fix_name(organism)
        name = organism.split()
        organism_name = " ".join(name[0:2])
        query = f'"{organism_name}"[Organism]'

        suffixes = name[2:]
        search_ids = entrez.read(
            entrez.esearch(db="genome", term=query, retmax=ENTREZ_RETMAX)
        )

        if not isinstance(search_ids, DictionaryElement):
            raise ValueError("Esearch should return a dictionary")

        search_results = entrez.read(
            entrez.esummary(db="genome", id=cast(Sequence[str], search_ids[K_ID_LIST]))
        )

        if search_results is None:
            return
        
        if not isinstance(search_results, ListElement):
            raise ValueError(f"Unexpected Entrez result: {search_results}")

        selected_organism: Optional[DictionaryElement] = None

        def choose_best(
            current: Optional[DictionaryElement],
            candidate: Optional[DictionaryElement]
        ):
            if current is None:
                return candidate
            
            if candidate is None:
                return current
            
            candidate_name = candidate[K_ORGANISM_NAME]
            assert isinstance(candidate_name, (str, StringElement)), "Unexpected organism name from entrez"

            for suffix in candidate_name[2:]:
                if suffix in suffixes:
                    return candidate

            current_name = current[K_ORGANISM_NAME]
            assert isinstance(current_name, (str, StringElement)), "Unexpected organism name from entrez"

            if len(candidate_name) < len(current_name):
                return candidate
            else:
                return current


        for item in search_results:
            if not isinstance(item, DictionaryElement):
                continue

            selected_organism = choose_best(selected_organism, item)

        if selected_organism is None:
            return

        accession = selected_organism[K_ACCESSION]

        assert isinstance(accession, (str, StringElement)), "Expected accession to be a string."

        for i,sequence in enumerate(self.__download_accession_gnomes(accession)):

            if organism_id is not None:
                sequence.id = f"{organism_id}_{i}"

            # todo: add the taxid
            yield SeqEntry(
                seq = sequence
            )
    
    def __load_from_sheet(self, file: str, wb: Workbook, spec: ExcelSheetSpec) -> Iterable[SeqEntryResult]:

        sheet = next(
            (wb[sheet] for sheet in wb.sheetnames if fnmatch(sheet, spec.name_pattern)),
            None
        )

        if sheet is None:
            return
        
        organism_records: Sequence[Cell] = sheet[spec.organism_column]
        id_records: Iterable[Optional[Cell]] = sheet[spec.id_column] \
                if spec.id_column is not None \
                else [None for _ in range(len(organism_records))]
        
        for orgn_cell, id_cell in zip(organism_records, id_records):
            orgn_name = orgn_cell.value
            orgn_id = id_cell.value if id_cell is not None else None

            if orgn_name is not None:
                try:
                    yield from self.__load_from_organism(
                        str(orgn_name),
                        str(orgn_id) if orgn_id is not None else None
                    )
                except Exception as e:
                    yield SequenceLoadException(
                        file=file,
                        message=f"Failed to load genome for \"{orgn_name}\": {e}"
                    )
        
    def __load_from_file(self, file_path: str, spec: ExcelSpec) -> Iterable[SeqEntryResult]:

        wb: Optional[Workbook] = None
        try:
            wb = openpyxl.open(file_path) #, read_only=True)
            for sheet_spec in spec.sheets:
                yield from self.__load_from_sheet(file_path, wb, sheet_spec)
        finally:
            wb.close() if wb is not None else None

class ExcelLoaderFactory(SeqLoaderFactoryBase):

    def __load_excel_specs(self, config_context: ConfigContext) -> Iterable[ExcelSpecConfig]:

        for config in config_context.configurations:
            try:
                spec = ExcelSpec(**config.configuration)
                yield ExcelSpecConfig(
                    spec = spec,
                    directory = path.dirname(config.path)
                )
            except ValidationError:
                # todo: advise if the type is excel
                pass

    def configure(self, config: ConfigContext) -> SeqLoaderBase:
        return ExcelLoader(
            list(self.__load_excel_specs(config))
        )
class Exporter:

    E_JSON = "JSON"
    E_BED = "BED"
    E_GFF3 = "GFF3"
    E_CSV = "CSV"

    def __new__(cls, filename: str, format: int = 0) -> None:
        if cls == Exporter:
            if format == Exporter.E_JSON:
                instance = super().__new__(JSON_Exporter)
            elif format == Exporter.E_BED:
                instance = super().__new__(BED_Exporter)
            elif format == Exporter.E_GFF3:
                instance = super().__new__(GFF3_Exporter)
            elif format == Exporter.E_CSV:
                instance = super().__new__(CSV_Exporter)
            else:
                print("ERROR: Invalid export format.")
        else:
            instance = super().__new__(cls)
        return instance
    
    def __init__(self, filename: str, format: int = 0):
        self.filename = filename

    
from .json_exporter import JSON_Exporter
from .bed_exporter import BED_Exporter
from .gff3_exporter import GFF3_Exporter
from .csv_exporter import CSV_Exporter

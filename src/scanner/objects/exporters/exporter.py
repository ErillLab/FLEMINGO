class Exporter:

    E_JSON = "JSON"
    E_BED = "BED"
    E_GFF3 = "GFF3"
    E_CSV = "CSV"

    def __new__(cls, filename, format=0):
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
    
    def __init__(self, filename, format=0):
        self.filename = filename

    
from FMS.objects.exporters.json_exporter import JSON_Exporter
from FMS.objects.exporters.bed_exporter import BED_Exporter
from FMS.objects.exporters.gff3_exporter import GFF3_Exporter
from FMS.objects.exporters.csv_exporter import CSV_Exporter

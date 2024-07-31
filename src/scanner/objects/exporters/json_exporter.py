from FMS.objects.exporters.exporter import Exporter

class JSON_Exporter(Exporter):

    def export(self, genome, placements_forw, placements_back):
        print("JSON exporter")


class FastaFile:
    def __init__(self,filename):
        self.filename = filename

    def __iter__(self):
        with open(self.filename,'r') as fh:
            for line in fh:
                line_s = line.strip()
                parse = self.parse_line(line_s)
                if parse:
                    yield parse
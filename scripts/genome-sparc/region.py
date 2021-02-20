import argparse
from collections import namedtuple
import re

# A named tuple that holds a chromosome:start-end triple:
class GRegion(namedtuple('GenomeRegion', ['chromosome', 'start', 'end'])):
    str_regex = re.compile('^\s*(?P<chromosome>[^:]+)\s*:\s*(?P<start>\d+)(?:\s*-\s*(?P<end>\d+))?\s*') # ^start of string \*s any number of white space ()group1 ?P call group chromosome
    @classmethod
    def fromString(cls, s):
        match = cls.str_regex.match(s)
        if match is None: raise ValueError('invalid chromosome region: {}'.format(s))
        if match.group('end') is None: return cls(match.group('chromosome'), int(match.group('start')), None)
        else: return cls(match.group('chromosome'), int(match.group('start')), int(match.group('end')))
    def flank(self, flank):
        return GRegion(self.chromosome, max(self.start - flank, 1), self.end + flank)
    def __str__(self):
        if self.end is None: return '{}:{}'.format(self.chromosome, self.start)
        return '{}:{}-{}'.format(self.chromosome, self.start, self.end)
    def __len__(self): return self.end - self.start


"""AGP v2.1 format reader and writer — adapted from RagTag ragtag_utilities/AGPFile.py.

Spec: https://www.ncbi.nlm.nih.gov/assembly/AGP_Specification/
"""
from dataclasses import dataclass


@dataclass
class AGPSeqLine:
    """AGP sequence component line (types: A/D/F/G/O/P/W)."""
    object_name: str
    object_beg: int
    object_end: int
    part_number: int
    component_type: str
    component_id: str
    component_beg: int
    component_end: int
    orientation: str  # '+', '-', '?', '0', 'na'

    def format(self) -> str:
        return "\t".join(str(x) for x in [
            self.object_name, self.object_beg, self.object_end,
            self.part_number, self.component_type, self.component_id,
            self.component_beg, self.component_end, self.orientation,
        ])


@dataclass
class AGPGapLine:
    """AGP gap line (types: N/U)."""
    object_name: str
    object_beg: int
    object_end: int
    part_number: int
    component_type: str
    gap_length: int
    gap_type: str
    linkage: str
    linkage_evidence: str

    def format(self) -> str:
        return "\t".join(str(x) for x in [
            self.object_name, self.object_beg, self.object_end,
            self.part_number, self.component_type, self.gap_length,
            self.gap_type, self.linkage, self.linkage_evidence,
        ])


class AGPReader:
    """Parse AGP v2.1 text."""

    VALID_SEQ_TYPES = {'A', 'D', 'F', 'G', 'O', 'P', 'W'}
    VALID_GAP_TYPES = {'N', 'U'}

    def parse(self, text: str) -> list:
        """Parse AGP text, returning AGPSeqLine and AGPGapLine objects."""
        lines = []
        for line in text.strip().split('\n'):
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) != 9:
                continue
            comp_type = fields[4]
            if comp_type in self.VALID_SEQ_TYPES:
                lines.append(AGPSeqLine(
                    object_name=fields[0],
                    object_beg=int(fields[1]),
                    object_end=int(fields[2]),
                    part_number=int(fields[3]),
                    component_type=comp_type,
                    component_id=fields[5],
                    component_beg=int(fields[6]),
                    component_end=int(fields[7]),
                    orientation=fields[8],
                ))
            elif comp_type in self.VALID_GAP_TYPES:
                lines.append(AGPGapLine(
                    object_name=fields[0],
                    object_beg=int(fields[1]),
                    object_end=int(fields[2]),
                    part_number=int(fields[3]),
                    component_type=comp_type,
                    gap_length=int(fields[5]),
                    gap_type=fields[6],
                    linkage=fields[7],
                    linkage_evidence=fields[8],
                ))
        return lines


class AGPWriter:
    """Write AGP v2.1 formatted output."""

    def format(self, lines: list) -> str:
        """Format a list of AGPLine objects to AGP text."""
        return '\n'.join(line.format() for line in lines) + '\n'

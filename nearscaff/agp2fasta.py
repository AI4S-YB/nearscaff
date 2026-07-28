"""Build scaffold FASTA from an AGP file plus the query contig FASTA.

Scaffolds described by the AGP are assembled first (reverse-complementing
components with '-' orientation and filling N/U gaps with Ns).  Query
contigs that were not placed in any scaffold are then appended verbatim
(full original header line preserved) at the end of the output file.
"""

_COMPLEMENT = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def collect_placed_components(agp_path: str) -> set:
    """Return the set of component ids placed in scaffolds.

    Component ids are column 6 of AGP lines whose type (column 5) is
    not a gap type (N/U).
    """
    placed = set()
    with open(agp_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            if fields[4] not in ('N', 'U'):
                placed.add(fields[5])
    return placed


def _read_agp_structure(agp_path: str) -> list:
    """Parse AGP into [(scaffold_name, [(kind, payload), ...]), ...].

    Each element is either ('gap', length) or
    ('seq', component_id, beg, end, orientation).
    """
    scaffolds = []
    current_name = None
    current_parts = []
    with open(agp_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9:
                continue
            if p[0] != current_name:
                if current_name is not None:
                    scaffolds.append((current_name, current_parts))
                current_name = p[0]
                current_parts = []
            if p[4] in ('N', 'U'):
                current_parts.append(('gap', int(p[5])))
            else:
                current_parts.append(
                    ('seq', p[5], int(p[6]), int(p[7]), p[8]))
    if current_name is not None:
        scaffolds.append((current_name, current_parts))
    return scaffolds


def _write_record(out, header: str, seq: str, line_width: int = 60):
    out.write(header if header.startswith('>') else '>' + header)
    out.write('\n')
    for i in range(0, len(seq), line_width):
        out.write(seq[i:i + line_width] + '\n')


def agp_to_fasta(agp_path: str, query_fasta: str, out_fasta: str,
                 gap_char: str = 'N', line_width: int = 60) -> tuple:
    """Write scaffold FASTA from AGP + query contigs.

    Scaffolds are written in AGP order.  Query contigs not present in the
    AGP are appended afterwards with their original (full) header lines.

    The query FASTA is streamed twice so that unplaced contigs are never
    held in memory; only sequences of placed components are loaded.

    Returns (n_scaffolds, n_unplaced_contigs).
    """
    placed = collect_placed_components(agp_path)
    scaffolds = _read_agp_structure(agp_path)

    # Pass 1: load sequences of placed components only
    contigs = {}
    name = None
    parts = []
    keep = False
    with open(query_fasta) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None and keep:
                    contigs[name] = ''.join(parts)
                name = line[1:].split()[0]
                keep = name in placed
                parts = []
            elif keep:
                parts.append(line.strip())
        if name is not None and keep:
            contigs[name] = ''.join(parts)

    with open(out_fasta, 'w') as out:
        # Scaffolds
        for scaf_name, parts_list in scaffolds:
            chunks = []
            for part in parts_list:
                if part[0] == 'gap':
                    chunks.append(gap_char * part[1])
                else:
                    _, cid, cbeg, cend, orient = part
                    seq = contigs.get(cid)
                    if seq is None:
                        # Component missing from query FASTA: fill with Ns
                        seq = gap_char * (cend - cbeg + 1)
                    else:
                        seq = seq[cbeg - 1:cend]
                        if orient == '-':
                            seq = _revcomp(seq)
                    chunks.append(seq)
            _write_record(out, scaf_name, ''.join(chunks), line_width)

        # Pass 2: append unplaced contigs verbatim (full header preserved)
        n_unplaced = 0
        write = False
        with open(query_fasta) as f:
            for line in f:
                if line.startswith('>'):
                    write = line[1:].split()[0] not in placed
                    if write:
                        n_unplaced += 1
                if write:
                    out.write(line if line.endswith('\n') else line + '\n')

    return len(scaffolds), n_unplaced

"""Unit tests for short-read PE-span gap size estimation."""

import pytest

from nearscaff.gapfill import parse_pe_spans


def _pe_paf(tmp_path, name, pairs):
    """pairs: [(read, te1_on_L, ts2_on_R)]; flank len 1500, reads 150 bp."""
    lines = []
    for rn, te1, ts2 in pairs:
        # R1 '+' on g1|L ending at te1 (150 bp aln)
        lines.append(f"{rn}\t150\t0\t150\t+\tg1|L\t1500\t{te1 - 150}\t{te1}"
                     f"\t150\t150\t60\ttp:A:P")
        # R2 '-' on g1|R starting at ts2
        lines.append(f"{rn}\t150\t0\t150\t-\tg1|R\t1500\t{ts2}\t{ts2 + 150}"
                     f"\t150\t150\t60\ttp:A:P")
    p1 = tmp_path / name
    p1.write_text("\n".join(l for i, l in enumerate(lines) if i % 2 == 0) + "\n")
    p2 = tmp_path / (name + ".2")
    p2.write_text("\n".join(l for i, l in enumerate(lines) if i % 2 == 1) + "\n")
    return str(p1), str(p2)


def test_pe_span_estimates_size(tmp_path):
    # two consistent pairs: shortfall 50 each side, 150 bp alignments
    # portion = 150 + 50 = 200 per mate -> est = 500 - 400 = 100
    p1, p2 = _pe_paf(tmp_path, "a.paf",
                     [("r1", 1450, 50), ("r2", 1440, 55)])
    est = parse_pe_spans(p1, p2, insert=500)
    assert "g1" in est
    assert 95 <= est["g1"] <= 110


def test_pe_span_requires_two_pairs(tmp_path):
    p1, p2 = _pe_paf(tmp_path, "a.paf", [("r1", 1450, 50)])
    assert parse_pe_spans(p1, p2, insert=500) == {}


def test_pe_span_rejects_wrong_orientation(tmp_path):
    # both reads '+' -> not a proper FR span
    lines = [
        "r1\t150\t0\t150\t+\tg1|L\t1500\t1300\t1450\t150\t150\t60\ttp:A:P",
        "r1\t150\t0\t150\t+\tg1|R\t1500\t50\t200\t150\t150\t60\ttp:A:P",
        "r2\t150\t0\t150\t+\tg1|L\t1500\t1300\t1450\t150\t150\t60\ttp:A:P",
        "r2\t150\t0\t150\t+\tg1|R\t1500\t50\t200\t150\t150\t60\ttp:A:P",
    ]
    p1 = tmp_path / "a.paf"
    p1.write_text(lines[0] + "\n" + lines[2] + "\n")
    p2 = tmp_path / "a.paf.2"
    p2.write_text(lines[1] + "\n" + lines[3] + "\n")
    assert parse_pe_spans(str(p1), str(p2), insert=500) == {}


def test_pe_span_rejects_low_mapq_and_secondary(tmp_path):
    p1, p2 = _pe_paf(tmp_path, "a.paf",
                     [("r1", 1450, 50), ("r2", 1450, 50)])
    # strip the tp tag and drop mapq -> all rejected
    p1x = tmp_path / "b.paf"
    p1x.write_text(open(p1).read().replace("\ttp:A:P", "").replace("\t60\n",
                                                                   "\t5\n"))
    p2x = tmp_path / "b.paf.2"
    p2x.write_text(open(p2).read().replace("\ttp:A:P", "").replace("\t60\n",
                                                                   "\t5\n"))
    assert parse_pe_spans(str(p1x), str(p2x), insert=500) == {}


def test_parse_pe_spans_preparsed_records_matches_file_path(tmp_path):
    # the sr branch now pre-parses each PAF once and passes records= ; confirm
    # that path yields identical results to the file-path path.
    from nearscaff.gapfill import _read_flank_paf
    p1, p2 = _pe_paf(tmp_path, "a.paf", [("r1", 1450, 50), ("r2", 1440, 55)])
    via_file = parse_pe_spans(p1, p2, insert=500)
    via_recs = parse_pe_spans(p1, p2, insert=500,
                              records1=_read_flank_paf(p1),
                              records2=_read_flank_paf(p2))
    assert via_file == via_recs
    assert "g1" in via_recs


# --- sr flow scaling: fetch only usable reads, scan each FASTQ once ---
# (root cause of the 10x hang: the old sr path fetched sequences for EVERY
#  read hitting any flank and re-scanned each FASTQ 4x via slow Python gzip.)

def _sr_paf_line(qn, ql, qs, qe, strand, tn, tl, ts, te):
    return (f"{qn}\t{ql}\t{qs}\t{qe}\t{strand}\t{tn}\t{tl}\t{ts}\t{te}"
            f"\t150\t150\t60\ttp:A:P")


def _sr_paf(tmp_path, name, lines):
    p = tmp_path / name
    p.write_text("".join(l + "\n" for l in lines))
    return str(p)


def _sr_fq(tmp_path, name, recs):
    out = []
    for rn, seq in recs.items():
        out.append(f"@{rn}\n{seq}\n+\n{'I' * len(seq)}\n")
    p = tmp_path / name
    p.write_text("".join(out))
    return str(p)


def test_sr_fill_fetches_only_usable_reads(tmp_path, monkeypatch):
    from nearscaff.gapfill import _sr_fill
    flank_len = {"g1|L": 1500, "g1|R": 1500}
    paf1 = _sr_paf(tmp_path, "r1.paf", [
        _sr_paf_line("rLonly", 150, 0, 150, "+", "g1|L", 1500, 1300, 1450),
        _sr_paf_line("rSpan", 150, 0, 150, "+", "g1|L", 1500, 1300, 1450),
        _sr_paf_line("rSpan", 150, 0, 150, "+", "g1|R", 1500, 50, 200),
        _sr_paf_line("rTail", 500, 0, 150, "+", "g1|L", 1500, 1385, 1495),
    ])
    paf2 = _sr_paf(tmp_path, "r2.paf", [])
    fq1 = _sr_fq(tmp_path, "r1.fq",
                 {"rLonly": "A" * 150, "rSpan": "A" * 150, "rTail": "A" * 500})
    fq2 = _sr_fq(tmp_path, "r2.fq", {})

    fetched = []

    def fake_fetch(fastqs, wanted):
        key = fastqs[0] if isinstance(fastqs, list) else fastqs
        fetched.append((key, set(wanted)))
        return {rn: "A" * 300 for rn in wanted}

    monkeypatch.setattr("nearscaff.gapfill._fetch_read_seqs", fake_fetch)
    monkeypatch.setattr("nearscaff.gapfill.overlap_closure_fills",
                        lambda *a, **k: {})

    _sr_fill(paf1, paf2, [fq1, fq2], flank_len, {"g1": "C" * 3000},
             str(tmp_path / "work"), threads=4)

    all_wanted = set().union(*(w for _, w in fetched)) if fetched else set()
    assert "rLonly" not in all_wanted, (
        f"single-flank non-tail read should not be fetched; got {all_wanted}")
    assert "rSpan" in all_wanted and "rTail" in all_wanted, (
        f"spanning + tail reads must be fetched; got {all_wanted}")


def test_sr_fill_scans_each_fastq_at_most_once(tmp_path, monkeypatch):
    from nearscaff.gapfill import _sr_fill
    flank_len = {"g1|L": 1500, "g1|R": 1500}
    paf1 = _sr_paf(tmp_path, "r1.paf", [
        _sr_paf_line("rSpan", 150, 0, 150, "+", "g1|L", 1500, 1300, 1450),
        _sr_paf_line("rSpan", 150, 0, 150, "+", "g1|R", 1500, 50, 200),
    ])
    paf2 = _sr_paf(tmp_path, "r2.paf", [
        _sr_paf_line("rSpan2", 150, 0, 150, "+", "g1|L", 1500, 1300, 1450),
        _sr_paf_line("rSpan2", 150, 0, 150, "+", "g1|R", 1500, 50, 200),
    ])
    fq1 = _sr_fq(tmp_path, "r1.fq", {"rSpan": "A" * 150})
    fq2 = _sr_fq(tmp_path, "r2.fq", {"rSpan2": "A" * 150})

    calls = []

    def fake_fetch(fastqs, wanted):
        calls.append(fastqs[0] if isinstance(fastqs, list) else fastqs)
        return {rn: "A" * 300 for rn in wanted}

    monkeypatch.setattr("nearscaff.gapfill._fetch_read_seqs", fake_fetch)
    monkeypatch.setattr("nearscaff.gapfill.overlap_closure_fills",
                        lambda *a, **k: {})

    _sr_fill(paf1, paf2, [fq1, fq2], flank_len, {"g1": "C" * 3000},
             str(tmp_path / "work"), threads=4)

    assert len(calls) <= 2, f"expected <=2 FASTQ scans, got {len(calls)}: {calls}"
    assert len(set(calls)) == len(calls), (
        f"same FASTQ scanned more than once: {calls}")


def test_sr_fill_produces_span_closure(tmp_path, monkeypatch):
    # regression guard: a read spanning both flanks of a small gap must yield
    # a span fill (s[b0:b1], here 50 bp). Must stay green through the stream
    # refactor -- behavior unchanged.
    from nearscaff.gapfill import _sr_fill
    flank_len = {"g1|L": 1500, "g1|R": 1500}
    paf1 = _sr_paf(tmp_path, "r1.paf", [
        # rSpan (500 bp) hits g1|L (qe=100, te=1450) and g1|R (qs=250, ts=50)
        # -> b0 = 100 + (1500-1450) = 150, b1 = 250 - 50 = 200 -> s[150:200] = 50bp
        _sr_paf_line("rSpan", 500, 0, 100, "+", "g1|L", 1500, 1350, 1450),
        _sr_paf_line("rSpan", 500, 250, 350, "+", "g1|R", 1500, 50, 150),
    ])
    paf2 = _sr_paf(tmp_path, "r2.paf", [])
    fq1 = _sr_fq(tmp_path, "r1.fq", {"rSpan": "ACGT" * 125})  # 500 bp, no N
    fq2 = _sr_fq(tmp_path, "r2.fq", {})
    monkeypatch.setattr("nearscaff.gapfill._fetch_read_seqs",
                        lambda fastqs, wanted: {rn: "ACGT" * 125 for rn in wanted})
    monkeypatch.setattr("nearscaff.gapfill.overlap_closure_fills",
                        lambda *a, **k: {})
    span, ovl, pe = _sr_fill(paf1, paf2, [fq1, fq2], flank_len,
                             {"g1": "C" * 3000}, str(tmp_path / "work"),
                             threads=4)
    assert "g1" in span, f"expected span closure for g1, got span={list(span)}"
    assert len(span["g1"]) == 50


def test_sr_fill_does_not_materialize_paf_records(tmp_path, monkeypatch):
    # stream contract: _sr_fill must NOT load the whole flank PAF into a
    # resident records list (the ~86GB bomb on 325M-read inputs). It should
    # stream each PAF once, building span/PE/tail/usable structures inline.
    from nearscaff import gapfill
    flank_len = {"g1|L": 1500, "g1|R": 1500}
    paf1 = _sr_paf(tmp_path, "r1.paf", [
        _sr_paf_line("rSpan", 150, 0, 150, "+", "g1|L", 1500, 1300, 1450),
        _sr_paf_line("rSpan", 150, 0, 150, "+", "g1|R", 1500, 50, 200),
    ])
    paf2 = _sr_paf(tmp_path, "r2.paf", [])
    fq1 = _sr_fq(tmp_path, "r1.fq", {"rSpan": "A" * 150})
    fq2 = _sr_fq(tmp_path, "r2.fq", {})

    materialized = []
    real = gapfill._read_flank_paf

    def spy(path):
        materialized.append(path)
        return real(path)

    monkeypatch.setattr(gapfill, "_read_flank_paf", spy)
    monkeypatch.setattr(gapfill, "_fetch_read_seqs",
                        lambda fastqs, wanted: {rn: "A" * 300 for rn in wanted})
    monkeypatch.setattr(gapfill, "overlap_closure_fills", lambda *a, **k: {})

    gapfill._sr_fill(paf1, paf2, [fq1, fq2], flank_len, {"g1": "C" * 3000},
                     str(tmp_path / "work"), threads=4)
    assert not materialized, (
        f"_sr_fill materialized the PAF into a records list via "
        f"_read_flank_paf ({materialized}); it must stream instead")


def _pe_pair_pafs(tmp_path):
    """Two primary FR pairs (rp1, rp2): R1 '+' on g1|L (te=1450), R2 '-' on
    g1|R (ts=50) -> with PE enabled yields pe_est['g1'] (~100 bp)."""
    flank_len = {"g1|L": 1500, "g1|R": 1500}
    paf1 = _sr_paf(tmp_path, "r1.paf", [
        _sr_paf_line("rp1", 150, 0, 150, "+", "g1|L", 1500, 1300, 1450),
        _sr_paf_line("rp2", 150, 0, 150, "+", "g1|L", 1500, 1300, 1450),
    ])
    paf2 = _sr_paf(tmp_path, "r2.paf", [
        _sr_paf_line("rp1", 150, 0, 150, "-", "g1|R", 1500, 50, 200),
        _sr_paf_line("rp2", 150, 0, 150, "-", "g1|R", 1500, 50, 200),
    ])
    return paf1, paf2, flank_len


def test_sr_fill_skips_pe_by_default(tmp_path, monkeypatch):
    # PE resize estimates are the collect-time bottleneck on huge inputs
    # (a 24GB pe_hits dict). sr mode skips PE by default -- pe_est must be
    # empty even when the input has valid FR spanning pairs.
    from nearscaff.gapfill import _sr_fill
    paf1, paf2, flank_len = _pe_pair_pafs(tmp_path)
    fq1 = _sr_fq(tmp_path, "r1.fq", {"rp1": "A" * 150, "rp2": "A" * 150})
    fq2 = _sr_fq(tmp_path, "r2.fq", {})
    monkeypatch.setattr("nearscaff.gapfill._fetch_read_seqs",
                        lambda f, w: {rn: "A" * 300 for rn in w})
    monkeypatch.setattr("nearscaff.gapfill.overlap_closure_fills",
                        lambda *a, **k: {})
    span, ovl, pe = _sr_fill(paf1, paf2, [fq1, fq2], flank_len,
                             {"g1": "C" * 3000}, str(tmp_path / "work"),
                             threads=4)
    assert pe == {}, f"sr must skip PE by default; got pe_est={pe}"


def test_sr_fill_with_pe_enabled_keeps_estimates(tmp_path, monkeypatch):
    # regression: with_pe=True still produces PE estimates (the capability
    # is preserved behind the flag).
    from nearscaff.gapfill import _sr_fill
    paf1, paf2, flank_len = _pe_pair_pafs(tmp_path)
    fq1 = _sr_fq(tmp_path, "r1.fq", {"rp1": "A" * 150, "rp2": "A" * 150})
    fq2 = _sr_fq(tmp_path, "r2.fq", {})
    monkeypatch.setattr("nearscaff.gapfill._fetch_read_seqs",
                        lambda f, w: {rn: "A" * 300 for rn in w})
    monkeypatch.setattr("nearscaff.gapfill.overlap_closure_fills",
                        lambda *a, **k: {})
    span, ovl, pe = _sr_fill(paf1, paf2, [fq1, fq2], flank_len,
                             {"g1": "C" * 3000}, str(tmp_path / "work"),
                             threads=4, with_pe=True)
    assert "g1" in pe, f"with_pe=True must yield PE estimate; got {pe}"


def test_sr_fill_caps_tail_reads_per_gap_side(tmp_path, monkeypatch):
    # The fetch bottleneck on huge inputs is usable being dominated by tail
    # reads (4.7M of 4.7M on the 10x run). Tails are capped at max_tails per
    # (gap, side) when built -- so collect must cap UPFRONT: only the first
    # max_tails tail reads per (gap, side) enter usable, the rest are never
    # fetched. Here 60 L-side tails on g1, max_tails=50 -> fetch <= 50.
    from nearscaff.gapfill import _sr_fill
    flank_len = {"g1|L": 1500, "g1|R": 1500}
    # 60 distinct tail reads on g1|L (ql=500, qe=150 -> clip 350 >= min_clip;
    # te=1495 -> within LR_EDGE=10 of the L edge)
    lines = [_sr_paf_line(f"t{i}", 500, 0, 150, "+", "g1|L", 1500, 1345, 1495)
             for i in range(60)]
    paf1 = _sr_paf(tmp_path, "r1.paf", lines)
    paf2 = _sr_paf(tmp_path, "r2.paf", [])
    fq1 = _sr_fq(tmp_path, "r1.fq", {f"t{i}": "A" * 500 for i in range(60)})
    fq2 = _sr_fq(tmp_path, "r2.fq", {})

    fetched = set()

    def fake_fetch(fastqs, wanted):
        fetched.update(wanted)
        return {rn: "A" * 500 for rn in wanted}

    monkeypatch.setattr("nearscaff.gapfill._fetch_read_seqs", fake_fetch)
    monkeypatch.setattr("nearscaff.gapfill.overlap_closure_fills",
                        lambda *a, **k: {})
    _sr_fill(paf1, paf2, [fq1, fq2], flank_len, {"g1": "C" * 3000},
             str(tmp_path / "work"), threads=4, max_tails=50)
    assert len(fetched) <= 50, (
        f"tail reads must be capped at max_tails=50 per (gap,side) BEFORE "
        f"fetch; fetched {len(fetched)}")


def test_fetch_read_seqs_handles_gz(tmp_path):
    # .gz inputs go through an external gzip -dc process (was Python gzip);
    # must still return the right {name: seq} for the wanted subset.
    import gzip as _gzip
    from nearscaff.gapfill import _fetch_read_seqs
    fq = tmp_path / "r.fq.gz"
    with _gzip.open(fq, "wt") as f:
        f.write("@r1\nACGTAC\n+\nIIIIII\n"
                "@r2\nGGGGGG\n+\nIIIIII\n"
                "@r3\nCCCCCC\n+\nIIIIII\n")
    seqs = _fetch_read_seqs(str(fq), {"r1", "r3"})
    assert seqs == {"r1": "ACGTAC", "r3": "CCCCCC"}


def test_sr_fill_parallel_matches_serial(tmp_path):
    # parallel=True (ProcessPoolExecutor, two worker processes) must give the
    # same span + overlap as serial. Uses a real FASTQ (no fetch mock -- the
    # mock would live in the parent, invisible to workers).
    from nearscaff.gapfill import _sr_fill
    flank_len = {"g1|L": 1500, "g1|R": 1500}
    paf1 = _sr_paf(tmp_path, "r1.paf", [
        _sr_paf_line("rSpan", 500, 0, 100, "+", "g1|L", 1500, 1350, 1450),
        _sr_paf_line("rSpan", 500, 250, 350, "+", "g1|R", 1500, 50, 150),
    ])
    paf2 = _sr_paf(tmp_path, "r2.paf", [])
    fq1 = _sr_fq(tmp_path, "r1.fq", {"rSpan": "ACGT" * 125})
    fq2 = _sr_fq(tmp_path, "r2.fq", {})
    kw = dict(reads_fq=[fq1, fq2], flank_len=flank_len,
              scaf_seqs={"g1": "C" * 3000}, work_prefix=str(tmp_path / "w"),
              threads=4)
    s_serial, o_serial, _ = _sr_fill(paf1, paf2, **kw)
    s_par, o_par, _ = _sr_fill(paf1, paf2, **kw, parallel=True)
    assert s_serial == s_par, f"span differs: {s_serial} vs {s_par}"
    assert o_serial == o_par, f"ovl differs: {o_serial} vs {o_par}"


# --- awk prefilter: C-speed line filtering feeding the Python parsers ---
# (collect spent ~18min of pure-Python PAF line parsing on the 10x lenta
# input; awk drops ~70% of lines at C speed so Python only parses the edge
# hits.  flank_len is preloaded into an awk array because synthetic/unit
# PAFs may carry a tl that differs from flank_len, and the tail test keys
# off flank_len, not tl.)

def _which_awk():
    import shutil
    return shutil.which("mawk") or shutil.which("awk")


def _edge_paf(tmp_path, name="edge.paf"):
    """PAF mixing every line class the collect prefilter must distinguish."""
    lines = [
        # span L hit: te 1450 >= tl-300=1200 (kept: span)
        _sr_paf_line("rSpanL", 150, 0, 150, "+", "g1|L", 1500, 1300, 1450),
        # span R hit: ts 50 <= 300 (kept: span)
        _sr_paf_line("rSpanR", 150, 0, 150, "+", "g1|R", 1500, 50, 200),
        # tail-only via flank_len != tl: flank_len says 1000 so
        # |te-1000|=5 <= 10 -> tail, but the tl-based span test (te>=1200)
        # would drop it -> exercises the flank_len preload (kept: tail)
        _sr_paf_line("rTailFl", 500, 0, 150, "+", "g2|L", 1500, 855, 1005),
        # mid-flank junk: low mapq, no tp tag, mid positions -> useless in
        # every mode (dropped).  NOTE a mid-flank primary mapq-60 line would
        # be legitimately kept with_pe=True (it feeds pe_hits); rPe covers
        # that; _sr_paf_line defaults would make this line PE-useful.
        "rJunk\t150\t0\t150\t+\tg1|L\t1500\t50\t200\t150\t150\t4",
        # primary mapq-60 mid-flank: PE-only (kept just with_pe=True)
        _sr_paf_line("rPe", 150, 0, 150, "+", "g1|L", 1500, 200, 350),
        # primary but mapq 4 mid-flank (dropped even with_pe)
        "rPeLow\t150\t0\t150\t+\tg1|L\t1500\t200\t350\t150\t150\t4\ttp:A:P",
        # secondary mapq-60 mid-flank (dropped even with_pe)
        "rPeSec\t150\t0\t150\t+\tg1|L\t1500\t200\t350\t150\t150\t60\ttp:A:S",
        # too few fields (ignored by awk NF>=12 and Python len<12 alike)
        "short\tline",
    ]
    flank_len = {"g1|L": 1500, "g1|R": 1500, "g2|L": 1000}
    return _sr_paf(tmp_path, name, lines), flank_len


def test_paf_edge_lines_filters_to_edge_hits(tmp_path):
    awk = _which_awk()
    if not awk:
        pytest.skip("no awk on PATH")
    import nearscaff.gapfill as gapfill
    paf, flank_len = _edge_paf(tmp_path)
    with gapfill._paf_edge_lines(paf, flank_len, with_pe=False) as f:
        kept = {l.split("\t")[0] for l in f}
    assert kept == {"rSpanL", "rSpanR", "rTailFl"}
    with gapfill._paf_edge_lines(paf, flank_len, with_pe=True) as f:
        kept_pe = {l.split("\t")[0] for l in f}
    assert kept_pe == {"rSpanL", "rSpanR", "rTailFl", "rPe"}


def test_sr_collect_matches_with_and_without_awk(tmp_path, monkeypatch):
    awk = _which_awk()
    if not awk:
        pytest.skip("no awk on PATH")
    import nearscaff.gapfill as gapfill
    paf, flank_len = _edge_paf(tmp_path)
    monkeypatch.setattr(gapfill, "_find_awk", lambda: awk)
    via_awk = gapfill._sr_collect(paf, flank_len, 30, 50, with_pe=True)
    monkeypatch.setattr(gapfill, "_find_awk", lambda: None)
    via_py = gapfill._sr_collect(paf, flank_len, 30, 50, with_pe=True)
    assert via_awk == via_py
    # non-trivial results: rPe is a primary mapq-60 mid-flank hit -> PE only
    assert "rPe" in via_py[1]
    assert "rPeLow" not in via_py[1] and "rPeSec" not in via_py[1]
    # rTailFl is usable only through the flank_len-keyed tail test
    assert "rTailFl" in via_py[3]


def test_sr_collect_empty_flank_len_matches(tmp_path, monkeypatch):
    # empty flank_len means an empty awk array file: the FILENAME-based
    # two-file idiom must not misread the PAF itself as array input
    awk = _which_awk()
    if not awk:
        pytest.skip("no awk on PATH")
    import nearscaff.gapfill as gapfill
    paf, _ = _edge_paf(tmp_path)
    monkeypatch.setattr(gapfill, "_find_awk", lambda: awk)
    via_awk = gapfill._sr_collect(paf, {}, 30, 50, with_pe=False)
    monkeypatch.setattr(gapfill, "_find_awk", lambda: None)
    via_py = gapfill._sr_collect(paf, {}, 30, 50, with_pe=False)
    assert via_awk == via_py


def test_fetch_read_seqs_header_space_dup_and_gz(tmp_path):
    # header "name extras" -> first token; duplicate name -> LAST wins;
    # .gz input -> gzip -dc subprocess.  (An awk fetch stage was benchmarked
    # and removed: the .gz path is decompression-bound, so the scan stays in
    # Python.)
    import gzip as _gzip

    from nearscaff.gapfill import _fetch_read_seqs
    body = ("@r1 extra 1:N:0:1\nACGTAC\n+\nIIIIII\n"
            "@r2\nGGGGGG\n+\nIIIIII\n"
            "@rDup\nAAAAAA\n+\nIIIIII\n"
            "@rDup\nTTTTTT\n+\nIIIIII\n"
            "@r3\nCCCCCC\n+\nIIIIII\n")
    fq = tmp_path / "x.fq"
    fq.write_text(body)
    gz = tmp_path / "y.fq.gz"
    with _gzip.open(gz, "wt") as f:
        f.write("@rGz\nATATAT\n+\nIIIIII\n")
    seqs = _fetch_read_seqs([str(fq), str(gz)], {"r1", "rDup", "r3", "rGz"})
    assert seqs == {"r1": "ACGTAC", "rDup": "TTTTTT", "r3": "CCCCCC",
                    "rGz": "ATATAT"}


def test_fetch_read_seqs_stops_between_files(tmp_path):
    # all wanted found in file1 -> file2 never read: r1 keeps file1's seq
    from nearscaff.gapfill import _fetch_read_seqs
    f1 = _sr_fq(tmp_path, "a.fq", {"r1": "AAAAAA"})
    f2 = _sr_fq(tmp_path, "b.fq", {"r1": "GGGGGG"})
    assert _fetch_read_seqs([f1, f2], {"r1"}) == {"r1": "AAAAAA"}

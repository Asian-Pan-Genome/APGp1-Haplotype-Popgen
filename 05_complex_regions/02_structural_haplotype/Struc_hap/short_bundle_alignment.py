"""Small, deterministic global MSA for bundle records below 100 bp."""


SHORT_BUNDLE_MAFFT_MIN_BP = 100


def align_short_bundle_records(records):
    """Align short sequences to the longest row without invoking MAFFT.

    Pairwise global alignments share one reference coordinate system. Inserted
    bases at each reference boundary are padded to the widest observed slot.
    This keeps every original residue and yields equal-width MSA rows.
    """
    if not records or any(not sequence or "-" in sequence for _, sequence in records):
        raise ValueError("short-bundle alignment requires nonempty ungapped sequences")
    if len({name for name, _ in records}) != len(records):
        raise ValueError("short-bundle alignment requires unique record IDs")
    reference = max(records, key=lambda item: len(item[1]))[1]
    reference_bp = len(reference)
    projections = {}
    insertion_widths = [0] * (reference_bp + 1)
    for name, sequence in records:
        n, m = reference_bp, len(sequence)
        scores = [[0] * (m + 1) for _ in range(n + 1)]
        trace = [[0] * (m + 1) for _ in range(n + 1)]
        for i in range(1, n + 1):
            scores[i][0], trace[i][0] = -2 * i, 1
        for j in range(1, m + 1):
            scores[0][j], trace[0][j] = -2 * j, 2
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                options = (scores[i - 1][j - 1] + (2 if reference[i - 1] == sequence[j - 1] else -1),
                           scores[i - 1][j] - 2, scores[i][j - 1] - 2)
                direction = max(range(3), key=lambda index: options[index])
                scores[i][j], trace[i][j] = options[direction], direction
        paired = []
        i, j = n, m
        while i or j:
            direction = trace[i][j]
            if direction == 0:
                paired.append((reference[i - 1], sequence[j - 1]))
                i, j = i - 1, j - 1
            elif direction == 1:
                paired.append((reference[i - 1], "-"))
                i -= 1
            else:
                paired.append(("-", sequence[j - 1]))
                j -= 1
        insertions = [""] * (reference_bp + 1)
        bases = []
        position = 0
        for reference_base, query_base in reversed(paired):
            if reference_base == "-":
                insertions[position] += query_base
            else:
                bases.append(query_base)
                position += 1
        if position != reference_bp:
            raise ValueError(f"invalid short-bundle reference projection: {name}")
        for boundary, insertion in enumerate(insertions):
            insertion_widths[boundary] = max(insertion_widths[boundary], len(insertion))
        projections[name] = insertions, bases
    aligned = {}
    for name, sequence in records:
        insertions, bases = projections[name]
        row = "".join(
            insertions[boundary].ljust(insertion_widths[boundary], "-")
            + (bases[boundary] if boundary < reference_bp else "")
            for boundary in range(reference_bp + 1)
        )
        if row.replace("-", "") != sequence:
            raise ValueError(f"short-bundle alignment changed residues: {name}")
        aligned[name] = row
    return aligned

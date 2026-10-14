VALID_BASES = {"A", "C", "G", "T"}


def parse_labeling_bases(labeling_base):
    bases = [b.strip().upper() for b in str(labeling_base).split(",") if b.strip()]
    if not bases:
        raise ValueError("--labeling_base must specify at least one base")
    bad = [b for b in bases if b not in VALID_BASES]
    if bad:
        raise ValueError(
            f"--labeling_base contains invalid base(s): {','.join(bad)}. Allowed: A,C,G,T"
        )
    if len(set(bases)) != len(bases):
        raise ValueError(f"--labeling_base contains duplicate base(s): {labeling_base}")
    return bases


def parse_conversion_map(labeling_base, sub_base):
    from_bases = parse_labeling_bases(labeling_base)

    if sub_base is None or str(sub_base).strip() == "":
        return from_bases, {}

    to_bases = [b.strip().upper() for b in str(sub_base).split(",") if b.strip()]
    bad = [b for b in to_bases if b not in VALID_BASES]
    if bad:
        raise ValueError(
            f"--sub_base contains invalid base(s): {','.join(bad)}. Allowed: A,C,G,T"
        )

    if len(to_bases) == 1 and len(from_bases) > 1:
        to_bases = to_bases * len(from_bases)

    if len(to_bases) != len(from_bases):
        raise ValueError(
            "--sub_base must have a single entry or exactly match the number of "
            f"--labeling_base entries ({len(from_bases)}); got {len(to_bases)}"
        )

    conv_map = {}
    for f, t in zip(from_bases, to_bases):
        if f == t:
            raise ValueError(
                f"--sub_base entry for labeling base '{f}' equals '{t}'; "
                "conversion must change the base"
            )
        conv_map[f] = t
    return from_bases, conv_map

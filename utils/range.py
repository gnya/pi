def split_range(start: int, end: int, parts: int) -> list[tuple[int, int]]:
    return [
        (
            start + ((end - start) * n) // parts,
            start + ((end - start) * (n + 1)) // parts,
        )
        for n in range(parts)
    ]

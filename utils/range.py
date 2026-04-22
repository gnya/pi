def split_range(size: int, parts: int) -> list[tuple[int, int]]:
    return [
        (
            (size * n) // parts,
            (size * (n + 1)) // parts - 1,
        )
        for n in range(parts)
    ]

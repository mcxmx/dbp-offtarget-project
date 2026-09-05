from __future__ import annotations

from collections import defaultdict
from typing import Iterable, Mapping


def stable_descending_order(values: Iterable[float]) -> list[int]:
    """Return deterministic descending order, preserving original order on ties."""
    indexed = list(enumerate(values))
    return [index for index, _ in sorted(indexed, key=lambda item: (-item[1], item[0]))]


def rank_order_unchanged_by_constant(values: Iterable[float], constant: float) -> bool:
    """Check the algebraic fact that subtracting a candidate-independent constant preserves ranking."""
    original = list(values)
    shifted = [value - constant for value in original]
    return stable_descending_order(original) == stable_descending_order(shifted)


class DisjointSet:
    def __init__(self, items: Iterable[str]) -> None:
        self.parent = {item: item for item in items}

    def find(self, item: str) -> str:
        parent = self.parent[item]
        if parent != item:
            self.parent[item] = self.find(parent)
        return self.parent[item]

    def union(self, left: str, right: str) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root != right_root:
            self.parent[right_root] = left_root

    def components(self) -> dict[str, list[str]]:
        grouped: dict[str, list[str]] = defaultdict(list)
        for item in sorted(self.parent):
            grouped[self.find(item)].append(item)
        return {root: members for root, members in grouped.items()}


def build_components(
    proteins: Iterable[str],
    pairwise_edges: Iterable[tuple[str, str]],
) -> dict[str, list[str]]:
    graph = DisjointSet(proteins)
    for left, right in pairwise_edges:
        graph.union(left, right)
    return graph.components()

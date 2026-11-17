"""The layer strings of one connected component, from its canonical form.

The layers are written on the canonical graph: the component relabeled by its
canonical numbers. Like the old AutoMol, the strings use 1-based numbers.
"""

from collections import Counter
from collections.abc import Sequence

import networkx as nx

from ._canon import CanonicalForm
from ._graph import MolGraph


def canonical_graph(graph: MolGraph, form: CanonicalForm) -> MolGraph:
    """The component of `form`, relabeled by its canonical numbers (0-based)."""
    component = graph.subgraph(form.numbers)
    return nx.relabel_nodes(component, form.numbers)


def formula_layer(graph: MolGraph) -> str:
    """The formula: C, then H, then the other symbols alphabetically; no count of 1."""
    counts = Counter(graph.nodes[k]["symbol"] for k in graph)
    counts["H"] += sum(graph.nodes[k]["hydrogens"] for k in graph)
    order = {"C": 0, "H": 1}
    symbols = sorted(
        (s for s, n in counts.items() if n), key=lambda s: (order.get(s, 2), s)
    )
    return "".join(s if counts[s] == 1 else f"{s}{counts[s]}" for s in symbols)


def c_layer(graph: nx.Graph) -> str:
    """The connection layer, without its prefix, of a canonical graph (D24).

    Written as InChI writes its connections layer (`GetDfsOrder4CT` in InChI's
    `ichimake.c`, called by `ichiprt3.c`):

    1. **Start atom:** the atom with the fewest neighbors, the lowest number on a tie.
    2. **Spanning tree:** a depth-first search from it that tries neighbors in
       ascending order. A neighbor already visited, other than the parent, that was
       visited earlier closes a ring; the closure belongs to the later atom.
    3. **Neighbor order:** each atom lists its ring closures first, then its children
       by the size of their subtrees, smallest first, where a ring closure inside the
       subtree counts as one more. Ties go to the lower number.
    4. **Writing:** the start atom, then each atom's listed neighbors in that order,
       a child followed at once by its own. A lone neighbor is preceded by `-`;
       otherwise the first by `(`, the last by `)` and the others by `,`. So the last
       child continues the chain: `1-2-3(4,5-6)7-8`, and cyclohexane is
       `1-2-4-6-5-3-1`.

    A graph with fewer than two atoms has no connection layer.
    """
    count = graph.number_of_nodes()
    if count <= 1:
        return ""
    neighbors = [sorted(graph[k]) for k in range(count)]
    start = min(range(count), key=lambda k: (len(neighbors[k]), k))

    # The depth-first search: the tree, the ring closures, and the subtree sizes
    visit_order = {start: 0}
    parents: dict[int, int] = {}
    children: list[list[int]] = [[] for _ in range(count)]
    closures: list[list[int]] = [[] for _ in range(count)]
    sizes = [1] * count
    stack = [(start, iter(neighbors[start]))]
    while stack:
        atom, remaining = stack[-1]
        for neighbor in remaining:
            if neighbor not in visit_order:
                visit_order[neighbor] = len(visit_order)
                parents[neighbor] = atom
                children[atom].append(neighbor)
                stack.append((neighbor, iter(neighbors[neighbor])))
                break
            if parents.get(atom) != neighbor and (
                visit_order[neighbor] < visit_order[atom]
            ):
                closures[atom].append(neighbor)
                sizes[atom] += 1
        else:
            stack.pop()
            if atom in parents:
                sizes[parents[atom]] += sizes[atom]

    # Writing: ring closures, then children by subtree size, then by number
    listed = [
        closures[k] + sorted(children[k], key=lambda c: (sizes[c], c))
        for k in range(count)
    ]
    parts = [f"{start + 1}"]
    writing = [(start, iter(enumerate(listed[start])))]
    while writing:
        atom, remaining = writing[-1]
        for position, neighbor in remaining:
            parts.append(_delimiter(position, len(listed[atom])))
            parts.append(f"{neighbor + 1}")
            if parents.get(neighbor) == atom:
                writing.append((neighbor, iter(enumerate(listed[neighbor]))))
            break
        else:
            writing.pop()
    return "".join(parts)


def _delimiter(position: int, count: int) -> str:
    """The delimiter before the neighbor at `position` of `count` listed neighbors."""
    if count == 1:
        return "-"
    if position == 0:
        return "("
    return ")" if position == count - 1 else ","


def h_layer(graph: MolGraph) -> str:
    """The hydrogen layer, without its prefix, of a canonical graph.

    Atoms are grouped by their number of implicit hydrogens, fewest first, and each
    run of consecutive numbers is written as a range: `2-3,8H,1,4-5H2`.
    """
    hydrogens = {k + 1: graph.nodes[k]["hydrogens"] for k in graph}
    groups = []
    for count in sorted(set(hydrogens.values()) - {0}):
        numbers = sorted(k for k, n in hydrogens.items() if n == count)
        ranges = ",".join(map(_range, _runs(numbers)))
        groups.append(f"{ranges}H" if count == 1 else f"{ranges}H{count}")
    return ",".join(groups)


def _runs(numbers: Sequence[int]) -> list[list[int]]:
    """Split sorted numbers into runs of consecutive ones."""
    runs: list[list[int]] = []
    for k in numbers:
        if runs and runs[-1][-1] == k - 1:
            runs[-1].append(k)
        else:
            runs.append([k])
    return runs


def _range(run: Sequence[int]) -> str:
    return f"{run[0]}-{run[-1]}" if len(run) > 1 else f"{run[0]}"

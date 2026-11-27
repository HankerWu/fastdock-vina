"""Exact finite-state sensitivity to contacts between variable Asn/Gln termini."""

from itertools import combinations

import numpy as np

from model import RT, pair_energy


def coupling(left, right):
    """Mixed finite difference of the four pair energies, with O then N coordinates."""
    distances = np.linalg.norm(left[:, None] - right[None, :], axis=-1)
    if np.min(distances) >= 8:
        return 0.0
    types = np.arange(2)  # deposited: O acceptor (0), N donor (1)
    radii = np.array([1.7, 1.8])

    def interaction(a, b):
        u, v = (types ^ a)[:, None], (types ^ b)[None, :]
        return pair_energy(distances, radii[u], radii[v], u != v).sum()

    return float(
        interaction(1, 1) - interaction(0, 1) - interaction(1, 0) + interaction(0, 0)
    )


class Component:
    """Compile a min-degree elimination order and cache the unliganded log partition."""

    def __init__(self, nodes, edges, field):
        self.nodes = sorted(nodes)
        self.field = np.asarray(field)[self.nodes]
        local = {v: i for i, v in enumerate(self.nodes)}
        self.pairs = [
            (
                tuple(sorted((local[i], local[j]))),
                np.array([[0.0, 0.0], [0.0, -jij / RT]]),
            )
            for (i, j), jij in edges.items()
            if i in local
        ]
        graph = {v: set() for v in local.values()}
        for (i, j), _ in self.pairs:
            graph[i].add(j)
            graph[j].add(i)
        self.order = []
        self.width = 0
        while graph:
            v = min(graph, key=lambda i: (len(graph[i]), i))
            neighbors = graph.pop(v)
            self.width = max(self.width, len(neighbors))
            for n in neighbors:
                graph[n].discard(v)
                graph[n].update(neighbors - {n})
            self.order.append(v)
        self.apo = self.log_partition(self.field)

    def log_partition(self, field):
        """Use log-space binary factors; the columns of every table follow sorted site IDs."""
        factors = self.pairs + [
            ((i,), np.array([0.0, -h / RT])) for i, h in enumerate(field)
        ]
        for v in self.order:
            adjacent = [f for f in factors if v in f[0]]
            factors = [f for f in factors if v not in f[0]]
            axes = sorted({i for dims, _ in adjacent for i in dims})
            total = sum(
                values.reshape([2 if i in dims else 1 for i in axes])
                for dims, values in adjacent
            )
            reduced = np.logaddexp.reduce(total, axis=axes.index(v))
            factors.append((tuple(i for i in axes if i != v), reduced))
        return float(sum(value for _, value in factors))


class CoupledAmides:
    """Binary energy sum h_i s_i + sum J_ij s_i s_j relative to the deposited labels."""

    def __init__(self, receptor):
        self.index = {site.residue: i for i, site in enumerate(receptor.sites)}
        self.edges = {}
        graph = {i: set() for i in self.index.values()}
        for i, j in combinations(graph, 2):
            jij = coupling(receptor.endpoints[i], receptor.endpoints[j])
            if jij != 0:
                self.edges[i, j] = jij
                graph[i].add(j)
                graph[j].add(i)
        self.components = []
        unseen = set(graph)
        while unseen:
            nodes = {min(unseen)}
            frontier = list(nodes)
            while frontier:
                new = graph[frontier.pop()] - nodes
                frontier.extend(new)
                nodes.update(new)
            unseen.difference_update(nodes)
            self.components.append(Component(nodes, self.edges, receptor.state_energy))

    def correction(self, details):
        change = np.zeros(len(self.index))
        for residue, _, ligand_delta, _ in details:
            change[self.index[residue]] = ligand_delta
        correction = 0.0
        for c in self.components:
            delta = change[c.nodes]
            if delta.any():
                correction += RT * (c.apo - c.log_partition(c.field + delta))
        return correction

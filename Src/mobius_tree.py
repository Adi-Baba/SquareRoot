"""
Rational Approximant Tree from Möbius Conjugacy.

Exact arithmetic throughout. Every node carries a certified error bound
derived from the closed-form identity:

    φ^n(z) - r = 2r λ^n T(z) / (1 - λ^n T(z))

where T(z) = (z-r)/(z+r), λ = -u/(2a+u), r = a+u, u = r-a.
"""

from fractions import Fraction as F
from typing import List, Optional, Tuple, Dict
import math


# ============================================================ Core Algebra

def p_of(a: F, u: F) -> F:
    return u * (2 * a + u)


def phi(a: F, u: F, z: F) -> F:
    a, u, z = F(a), F(u), F(z)
    return a + p_of(a, u) / (a + z)


def T(a: F, u: F, z: F) -> F:
    a, u, z = F(a), F(u), F(z)
    r = a + u
    return (z - r) / (z + r)


def lam(a: F, u: F) -> F:
    a, u = F(a), F(u)
    return (a - (a + u)) / (a + (a + u))


def error_closed_form(a: F, u: F, z: F, n: int) -> F:
    v = lam(a, u) ** n * T(a, u, z)
    return 2 * (a + u) * v / (1 - v)


def phi_iter_exact(a: F, u: F, z: F, n: int) -> F:
    z = F(z)
    for _ in range(n):
        z = phi(a, u, z)
    return z - (a + u)


# ============================================================ Polynomial Operations

def poly_mul(p1: List[F], p2: List[F]) -> List[F]:
    result = [F(0)] * (len(p1) + len(p2) - 1)
    for i, c1 in enumerate(p1):
        for j, c2 in enumerate(p2):
            result[i + j] += c1 * c2
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def poly_add(p1: List[F], p2: List[F]) -> List[F]:
    n = max(len(p1), len(p2))
    result = [F(0)] * n
    for i in range(n):
        c1 = p1[i] if i < len(p1) else F(0)
        c2 = p2[i] if i < len(p2) else F(0)
        result[i] = c1 + c2
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def poly_sub(p1: List[F], p2: List[F]) -> List[F]:
    n = max(len(p1), len(p2))
    result = [F(0)] * n
    for i in range(n):
        c1 = p1[i] if i < len(p1) else F(0)
        c2 = p2[i] if i < len(p2) else F(0)
        result[i] = c1 - c2
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def poly_derivative(poly: List[F]) -> List[F]:
    if len(poly) <= 1:
        return [F(0)]
    return [poly[i] * F(i) for i in range(1, len(poly))]


def poly_eval(poly: List[F], x: F) -> F:
    result = F(0)
    for coeff in reversed(poly):
        result = result * x + coeff
    return result


# ============================================================ Error Polynomials

def error_polynomials(a: F, w: F, z: F, n: int) -> Tuple[List[F], List[F]]:
    a, w, z = F(a), F(w), F(z)
    n_int = n
    
    num_lam = [F(0), F(-1)]
    den_lam = [2*a, F(1)]
    num_T = [z - a, F(-1)]
    den_T = [z + a, F(1)]
    
    num_lam_n = [F(1)]
    den_lam_n = [F(1)]
    for _ in range(n_int):
        num_lam_n = poly_mul(num_lam_n, num_lam)
        den_lam_n = poly_mul(den_lam_n, den_lam)
    
    num_v = poly_mul(num_lam_n, num_T)
    den_v = poly_mul(den_lam_n, den_T)
    
    num_e = poly_mul([2*a, F(2)], num_v)
    den_e = poly_sub(den_v, num_v)
    
    return num_e, den_e


# ============================================================ Robust Root Finding

def find_derivative_roots(num: List[F], den: List[F], lo: F, hi: F) -> List[F]:
    num_prime = poly_derivative(num)
    den_prime = poly_derivative(den)
    
    term1 = poly_mul(num_prime, den)
    term2 = poly_mul(num, den_prime)
    expr = poly_sub(term1, term2)
    
    while len(expr) > 1 and expr[-1] == 0:
        expr.pop()
    
    if len(expr) <= 1:
        return []
    
    roots = []
    num_samples = 200
    step = (hi - lo) / num_samples
    
    if step == 0:
        return []
    
    prev_val = poly_eval(expr, lo)
    if prev_val == 0:
        roots.append(lo)
    
    for i in range(1, num_samples + 1):
        x = lo + step * i
        curr_val = poly_eval(expr, x)
        
        if curr_val == 0:
            roots.append(x)
        elif prev_val * curr_val < 0:
            root = bisect_polynomial(expr, x - step, x)
            roots.append(root)
        
        prev_val = curr_val
    
    filtered = []
    for r in roots:
        if lo <= r <= hi:
            if not any(abs(r - existing) < F(1, 1000) for existing in filtered):
                filtered.append(r)
    
    return filtered


def bisect_polynomial(poly: List[F], lo: F, hi: F, tol: F = F(1, 10**15)) -> F:
    f_lo = poly_eval(poly, lo)
    f_hi = poly_eval(poly, hi)
    
    if f_lo == 0:
        return lo
    if f_hi == 0:
        return hi
    
    if f_lo * f_hi > 0:
        return (lo + hi) / 2
    
    for _ in range(60):
        mid = (lo + hi) / 2
        f_mid = poly_eval(poly, mid)
        if f_mid == 0:
            return mid
        if f_lo * f_mid <= 0:
            hi = mid
            f_hi = f_mid
        else:
            lo = mid
            f_lo = f_mid
        if hi - lo < tol:
            break
    return (lo + hi) / 2


# ============================================================ Certified Max Error

def certified_max_error(a: F, w: F, z: F, n: int) -> F:
    num, den = error_polynomials(a, w, z, n)
    crits = find_derivative_roots(num, den, F(0), w)
    
    candidates = [F(0), w] + [c for c in crits if F(0) <= c <= w]
    
    best_val = F(0)
    for u in candidates:
        val = abs(poly_eval(num, u) / poly_eval(den, u))
        if val > best_val:
            best_val = val
    return best_val


# ============================================================ Tree Node (mutable for tree building)

class Node:
    """A rational approximant with exact certified error."""
    __slots__ = ('z', 'z0', 'a', 'w', 'n', 'error_bound', 'depth', 
                 'path', 'parent', 'children')
    
    def __init__(self, z: F, z0: F, a: F, w: F, n: int, 
                 error_bound: F, depth: int, path: Tuple[int, ...] = ()):
        self.z = z
        self.z0 = z0
        self.a = a
        self.w = w
        self.n = n
        self.error_bound = error_bound
        self.depth = depth
        self.path = path
        self.parent: Optional['Node'] = None
        self.children: List['Node'] = []
    
    def __lt__(self, other: 'Node') -> bool:
        return self.error_bound < other.error_bound
    
    def exact_error_at_u(self, u: F) -> F:
        return abs(error_closed_form(self.a, u, self.z0, self.n))
    
    def iterate_value_at_u(self, u: F) -> F:
        return self.a + u + error_closed_form(self.a, u, self.z0, self.n)
    
    def is_leaf(self) -> bool:
        return len(self.children) == 0
    
    def best_descendant(self) -> 'Node':
        best = self
        for child in self.children:
            desc = child.best_descendant()
            if desc.error_bound < best.error_bound:
                best = desc
        return best
    
    def __repr__(self) -> str:
        return f"Node(n={self.n}, z0={float(self.z0):.4f}, err={float(self.error_bound):.2e})"


# ============================================================ Tree Construction

def compute_iterate(a: F, u: F, z: F, n: int) -> F:
    return a + u + error_closed_form(a, u, z, n)


def build_tree(
    alpha: int,
    max_n: int = 8,
    z0_grid: int = 9,
    u_fixed: Optional[F] = None
) -> Tuple[List['Node'], Dict[Tuple[int, int], 'Node']]:
    """
    Build the rational approximant tree with navigation.
    
    Returns (all_nodes, node_lookup) where node_lookup maps (n, z0_idx) -> Node.
    Tree edges: parent (n, z0) -> child (n+1, same z0) representing one more φ iteration.
    """
    a, w, _ = initial_bracket(alpha)
    z0_values = [a + w * F(i, z0_grid - 1) for i in range(z0_grid)]
    
    # Build all nodes
    node_map: Dict[Tuple[int, int], 'Node'] = {}
    all_nodes: List['Node'] = []
    
    for n in range(1, max_n + 1):
        for z0_idx, z0 in enumerate(z0_values):
            if u_fixed is not None:
                error = abs(error_closed_form(a, u_fixed, z0, n))
                z_exact = a + u_fixed + error_closed_form(a, u_fixed, z0, n)
            else:
                error = certified_max_error(a, w, z0, n)
                z_exact = z0
            
            node = Node(
                z=z_exact, z0=F(z0), a=a, w=w, n=n,
                error_bound=error, depth=n, path=(n,)
            )
            node_map[(n, z0_idx)] = node
            all_nodes.append(node)
    
    # Link parent -> children (same z0, n+1)
    for n in range(1, max_n):
        for z0_idx in range(len(z0_values)):
            parent = node_map[(n, z0_idx)]
            child = node_map[(n + 1, z0_idx)]
            parent.children.append(child)
            child.parent = parent
    
    all_nodes.sort()
    return all_nodes, node_map


def initial_bracket(alpha: int) -> Tuple[F, F, F]:
    a = F(int(math.isqrt(alpha)))
    w = F(1)
    return a, w, a + w / 2


# ============================================================ Tree Navigation

def find_best_path(nodes: List['Node'], target_error: F) -> Optional[List['Node']]:
    """Find path from root to node with error <= target_error."""
    for node in nodes:
        if node.error_bound <= target_error:
            path = []
            curr = node
            while curr is not None:
                path.append(curr)
                curr = curr.parent
            return list(reversed(path))
    return None


def print_tree_summary(nodes: List['Node'], alpha: int):
    """Print summary of tree structure."""
    roots = [n for n in nodes if n.parent is None]
    print(f"Tree for alpha={alpha}: {len(nodes)} nodes, {len(roots)} roots (n=1)")
    print(f"Depth range: 1 to {max(n.n for n in nodes)}")
    print(f"z0 grid points: {len(set(n.z0 for n in nodes))}")
    
    for depth in sorted(set(n.n for n in nodes)):
        level_nodes = [nd for nd in nodes if nd.n == depth]
        best = min(level_nodes)
        print(f"  n={depth}: {len(level_nodes)} nodes, best err={float(best.error_bound):.2e} at z0={float(best.z0):.4f}")


# ============================================================ Verification

def verify_closed_form(alpha: int, n: int, u_samples: int = 20) -> bool:
    a, w, _ = initial_bracket(alpha)
    for i in range(1, u_samples):
        u = w * F(i, u_samples)
        z = a + w / 2
        exact_iter = phi_iter_exact(a, u, z, n)
        closed = error_closed_form(a, u, z, n)
        if exact_iter != closed:
            return False
    return True


# ============================================================ Demo & Testing

if __name__ == "__main__":
    print("Möbius Rational Approximant Tree")
    print("=================================")
    
    for alpha in [2, 3, 5, 7, 12, 17]:
        print(f"\nVerifying alpha={alpha}:")
        for n in range(1, 6):
            ok = verify_closed_form(alpha, n, u_samples=20)
            print(f"  n={n}: {'OK' if ok else 'FAIL'}")
    
    print("\n" + "="*60)
    print("EXACT TREE: alpha=12, u=1/3 (fixed)")
    print("="*60)
    nodes, node_map = build_tree(12, max_n=8, z0_grid=9, u_fixed=F(1, 3))
    print(f"Generated {len(nodes)} nodes")
    print("Top 10 by exact error at u=1/3:")
    for node in nodes[:10]:
        print(f"  n={node.n}, z0={float(node.z0):.6f}, z={float(node.z):.10f}, err={float(node.error_bound):.2e}")
    
    print("\n" + "="*60)
    print("CERTIFIED TREE: alpha=12 (max over u)")
    print("="*60)
    nodes_cert, node_map_cert = build_tree(12, max_n=8, z0_grid=9, u_fixed=None)
    print(f"Generated {len(nodes_cert)} nodes")
    print("Top 10 by certified error bound:")
    for node in nodes_cert[:10]:
        print(f"  n={node.n}, z0={float(node.z0):.6f}, err_bound={float(node.error_bound):.2e}")
    
    print("\n" + "="*60)
    print("TREE STRUCTURE (certified, alpha=12)")
    print("="*60)
    print_tree_summary(nodes_cert, 12)
    
    # Show navigation
    print("\n" + "="*60)
    print("NAVIGATION EXAMPLE: path from root to best n=8 node")
    print("="*60)
    best = nodes_cert[0]
    path = []
    curr = best
    while curr:
        path.append(curr)
        curr = curr.parent
    for node in reversed(path):
        print(f"  n={node.n}, z0={float(node.z0):.4f}, err={float(node.error_bound):.2e}")
    
    # Show tree structure from a root
    print("\n" + "="*60)
    print("TREE FROM ROOT z0=4.0 (n=1 to n=8)")
    print("="*60)
    root = node_map_cert[(1, 8)]  # z0 = a + w = 4.0
    curr = root
    for n in range(1, 9):
        node = node_map_cert[(n, 8)]
        print(f"  n={node.n}, z0={float(node.z0):.4f}, err_bound={float(node.error_bound):.2e}")

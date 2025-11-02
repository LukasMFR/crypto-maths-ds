# ================================================
#  Bézout, inverse modulaire et congruences (NumWorks)
#  - Décomposition en facteurs premiers
#  - PGCD / Bézout (Euclide étendu + remontée)
#  - Inverse modulaire (avec étapes)
#  - Équations modulaires (ax ≡ b [m])
#  - Tables d’addition/multiplication en Z et Z_n (+ Z_n et Z_n*)
#  - Théorème des restes chinois (CRT) pour systèmes copremiers
#  - Puissance mod m
#  Auteur : toi ;)
#  Version menu "one-shot" (pas de boucle, version allégée)
# ================================================

# ---------- Outils d'affichage ----------
def sep(title=None):
    # Affichage compact sur une seule ligne (compatible petit écran)
    if title:
        print("--- {} ---".format(title))
    else:
        print("-" * 38)

def show_divisions(divs):
    for (A, B, q, r) in divs:
        print("{} = {}*{} + {}".format(A, q, B, r))

def _expr_to_string(expr, R):
    terms = []
    for i in range(len(R)):
        c = expr.get(i, 0)
        if c != 0:
            terms.append("{}*{}".format(c, R[i]))
    if not terms:
        return "0"
    return " + ".join(terms)

# --- remplaçant rjust (compat MicroPython NumWorks) ---
def _rjust(val, w):
    s = str(val)
    L = len(s)
    if L < w:
        return " " * (w - L) + s
    return s

def _col_width(values):
    w = 1
    for v in values:
        l = len(str(v))
        if l > w:
            w = l
    return w + 1  # petite marge

def _print_table(header_vals, row_vals, cell_fn, title):
    sep(title)
    candidates = list(header_vals) + list(row_vals)
    for i in row_vals:
        for j in header_vals:
            candidates.append(cell_fn(i, j))
    w = _col_width(candidates)

    # en-tête
    first_cell = " " * w + "|"
    print(first_cell, end="")
    for j in header_vals:
        print(_rjust(j, w), end="")
    print()

    # séparation
    print("-" * (w + 1 + w * len(header_vals)))

    # lignes
    for i in row_vals:
        print(_rjust(i, w) + "|", end="")
        for j in header_vals:
            v = cell_fn(i, j)
            print(_rjust(v, w), end="")
        print()

# ---------- gcd silencieux (pour Z_n*) ----------
def gcd(a, b):
    a = abs(a); b = abs(b)
    while b:
        a, b = b, a % b
    return a

# ---------- Décomposition en facteurs premiers (méthode "échelle") ----------
def _format_factorization(fdict):
    # fdict: {prime: exponent} -> "2^4 * 3 * 7"
    parts = []
    for p in sorted(fdict.keys()):
        e = fdict[p]
        if e == 1:
            parts.append(str(p))
        else:
            parts.append("{}^{}".format(p, e))
    return " * ".join(parts) if parts else "1"

def prime_factors_ladder(n):
    if n == 0:
        print("0 : factorisation non définie (multiple de tous les entiers).")
        return {}
    sign = -1 if n < 0 else 1
    if sign < 0:
        print("Attention: n < 0 -> on factorise |n| et on garde le signe -1.")
    n = abs(n)
    if n == 1:
        print("1")
        return {}

    print(n)  # ligne de départ de l'échelle
    f = {}

    # facteur 2
    while n % 2 == 0:
        n //= 2
        print("{} | {}".format(n, 2))
        f[2] = f.get(2, 0) + 1

    # facteurs impairs
    p = 3
    while p * p <= n:
        while n % p == 0:
            n //= p
            print("{} | {}".format(n, p))
            f[p] = f.get(p, 0) + 1
        p += 2

    # reste premier > 1
    if n > 1:
        print("1 | {}".format(n))
        f[n] = f.get(n, 0) + 1

    return f if sign > 0 else ({-1:1} | f) if hasattr(dict, "__or__") else (dict([(-1,1)]) | f)

# ---------- Euclide étendu avec traçage + remontée ----------
def egcd_verbose(a, b, show=True, show_back=True):
    A0, B0 = a, b
    divs = []
    R = [a, b]
    Q = [None, None]

    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1

    while r != 0:
        q = old_r // r
        rem = old_r - q * r
        divs.append((old_r, r, q, rem))
        R.append(rem)
        Q.append(q)
        old_r, r = r, rem
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t

    g, x, y = old_r, old_s, old_t
    k = len(R) - 2

    if show:
        sep("Algorithme d'Euclide ({} , {})".format(A0, B0))
        show_divisions(divs)
        print("pgcd({}, {}) = {}".format(A0, B0, g))

    if show and show_back and k >= 2:
        sep("Remontée (combinaison linéaire)")
        print("{} = {} - {}*{}".format(R[k], R[k-2], Q[k], R[k-1]))
        expr = {k-2: 1, k-1: -Q[k]}
        for j in range(k-1, 1, -1):
            cj = expr.get(j, 0)
            if cj == 0:
                continue
            print("Remplacer {} par {} - {}*{}".format(R[j], R[j-2], Q[j], R[j-1]))
            expr.pop(j, None)
            expr[j-2] = expr.get(j-2, 0) + cj
            expr[j-1] = expr.get(j-1, 0) - cj * Q[j]
            print("=> {} = {}".format(R[k], _expr_to_string(expr, R)))
        x_back = expr.get(0, 0)
        y_back = expr.get(1, 0)
        print("Donc {} = {}*{} + {}*{}".format(g, x_back, A0, y_back, B0))

    if show:
        print("Coeffs de Bézout : x = {}, y = {}".format(x, y))
        print("Vérif : {}*{} + {}*{} = {}".format(A0, x, B0, y, A0*x + B0*y))
    return g, x, y

# ---------- Inverse modulaire ----------
def inv_mod(a, m, show=True):
    if m <= 0:
        if show:
            print("Le module doit être > 0")
        return False, None
    g, x, y = egcd_verbose(a, m, show=show, show_back=True)
    if g != 1:
        if show:
            sep("Inverse modulaire")
            print("gcd({}, {}) = {} ≠ 1 : pas d'inverse modulo {}.".format(a, m, g, m))
        return False, None
    inv = x % m
    if show:
        sep("Inverse mod {}".format(m))
        print("Inverse trouvé : {}^(-1) ≡ {}  [ {} ]".format(a, inv, m))
        print("Vérif : ({}*{}) % {} = {}".format(a, inv, m, (a*inv) % m))
    return True, inv

# ---------- Résolution a x ≡ b [m] ----------
def solve_congruence(a, b, m, show=True, list_rep=True):
    if m <= 0:
        if show:
            print("Le module doit être > 0")
        return False, None, None, None
    sep("Résolution de {} x ≡ {}  [ {} ]".format(a, b, m))
    d, xg, yg = egcd_verbose(a, m, show=True, show_back=False)
    if b % d != 0:
        print("Comme {} ne divise pas {}, aucune solution.".format(d, b))
        return False, None, None, d
    a1 = a // d; b1 = b // d; m1 = m // d
    print("On réduit : a'={}, b'={}, m'={} (d = {})".format(a1, b1, m1, d))
    print("Nouvelle équation : {} x ≡ {}  [ {} ]".format(a1, b1, m1))
    ok, inv = inv_mod(a1, m1, show=True)
    if not ok:
        print("Problème inattendu : a' et m' ne sont pas copremiers.")
        return False, None, None, d
    x0 = (inv * b1) % m1
    sep("Solution")
    print("Solution de base : x0 ≡ {} * {} ≡ {}  [ {} ]".format(inv, b1, x0, m1))
    print("Vérif : ({}*{}) % {} = {}  (doit ≡ {})".format(a, x0, m, (a*x0) % m, b % m))
    print("Forme générale des solutions : x ≡ {}  [ {} ]".format(x0, m1))
    if d > 1:
        k_values = ", ".join(str(i) for i in range(d))
        print("Remontée au mod {} : x ≡ {} + {}k, k = {}.".format(m, x0, m1, k_values))
        reps = []
        for k in range(d):
            reps.append((x0 + k*m1) % m)
        reps = sorted(list(set(reps)))
        if list_rep:
            print("Représentants (mod {}): {}".format(m, ", ".join(str(r) for r in reps)))
    return True, x0, m1, d

# ---------- Résolution d'équations : wrappers pour formes générales ----------

def solve_ax_plus_c_eq_b(a, c, b, m):
    sep("0) Mise en forme / réduction modulo {}".format(m))
    if (b % m) == 0 and (c % m) != 0:
        # Cas typique "a x + c ≡ 0"
        print("{} x ≡ -{} ≡ {}  [ {} ]".format(a, c, (-c) % m, m))
    elif (c % m) == 0:
        print("{} x ≡ {}  [ {} ]".format(a, b % m, m))
    else:
        print("{} x ≡ {} - {} ≡ {}  [ {} ]".format(a, b % m, c % m, (b - c) % m, m))
    return solve_congruence(a, (b - c) % m, m, show=True, list_rep=True)

def solve_x_plus_c_eq_b(c, b, m):
    sep("0) Mise en forme / réduction modulo {}".format(m))
    print("{} ≡ {}  [ {} ]".format(c, c % m, m))
    print("{} ≡ {}  [ {} ]".format(b, b % m, m))
    print("=> x ≡ {}  [ {} ]".format((b - c) % m, m))
    return solve_congruence(1, (b - c) % m, m, show=True, list_rep=True)

def equations_menu_option():
    sep("Équations - choisir une forme")
    print("1) a x ≡ b  [ m ]")
    print("2) a x + c ≡ 0  [ m ]")
    print("3) a x + c ≡ b  [ m ]")
    print("4) x + c ≡ b  [ m ]")
    ch = input("> Choix forme : ").strip()

    try:
        if ch == "1":
            a = int(input("a = "))
            b = int(input("b = "))
            m = int(input("m (module > 0) = "))
            solve_congruence(a, b, m, show=True, list_rep=True)

        elif ch == "2":
            a = int(input("a = "))
            c = int(input("c = "))
            m = int(input("m (module > 0) = "))
            sep("Équation {}x + {} = 0  [ {} ]".format(a, c, m))
            print("Réécriture : {}x ≡ -{} ≡ {}  [ {} ]".format(a, c, (-c) % m, m))
            solve_congruence(a, (-c) % m, m, show=True, list_rep=True)

        elif ch == "3":
            a = int(input("a = "))
            c = int(input("c = "))
            b = int(input("b = "))
            m = int(input("m (module > 0) = "))
            solve_ax_plus_c_eq_b(a, c, b, m)

        elif ch == "4":
            c = int(input("c = "))
            b = int(input("b = "))
            m = int(input("m (module > 0) = "))
            solve_x_plus_c_eq_b(c, b, m)

        else:
            print("Choix inconnu.")
    except:
        print("Entrée invalide.")

# ---------- Tables Z / Z_n ----------
def table_Z(start, end, op):
    if start > end:
        start, end = end, start
    rows = list(range(start, end + 1))
    cols = list(range(start, end + 1))
    if op == "+":
        cell = lambda i, j: i + j
        title = "Table d'addition en Z, [{}..{}]".format(start, end)
    else:
        cell = lambda i, j: i * j
        title = "Table de multiplication en Z, [{}..{}]".format(start, end)
    _print_table(cols, rows, cell, title)

def table_Zn(n, op):
    if n <= 0:
        print("Le module n doit être > 0")
        return
    rows = list(range(0, n))
    cols = list(range(0, n))
    if op == "+":
        cell = lambda i, j: (i + j) % n
        title = "Table d'addition modulo {}".format(n)
    else:
        cell = lambda i, j: (i * j) % n
        title = "Table de multiplication modulo {}".format(n)
    _print_table(cols, rows, cell, title)

    # Ensembles Z_n et Z_n*
    Zn_list = [i for i in range(n)]
    print("\nZ_{} = {{ {} }}".format(n, ", ".join([str(x) for x in Zn_list])))
    Zn_star = [a for a in range(n) if gcd(a, n) == 1]
    print("Z_{}* = {{ {} }}".format(n, ", ".join([str(x) for x in Zn_star])))

# ---------- CRT (méthode du prof - formule directe) ----------
def solve_system_crt_coprime():
    sep("CRT (formule directe)")
    k = int(input("Nombre d'équations k = "))
    if k <= 0:
        print("k doit être >= 1")
        return

    residues = []
    moduli = []
    for i in range(1, k+1):
        print("Equation #{} :".format(i))
        ai = int(input("  a{} = ".format(i)))
        mi = int(input("  m{} (>0) = ".format(i)))
        if mi <= 0:
            print("Module > 0 requis.")
            return
        ai = ai % mi
        residues.append(ai)
        moduli.append(mi)

    # Vérif coprimalité par paires
    for i in range(k):
        for j in range(i+1, k):
            if gcd(moduli[i], moduli[j]) != 1:
                print("Moduli NON copremiers (m{}={}, m{}={}).".format(i+1, moduli[i], j+1, moduli[j]))
                print("Utilise l'option (6) Système modulaire (cas général).")
                return

    # Système
    print("Système :")
    for i in range(k):
        print("x ≡ {}  [ {} ]".format(residues[i], moduli[i]))

    # Produit total
    M = 1
    for mi in moduli:
        M *= mi
    sep("Produit total")
    prod_str = " * ".join(str(mi) for mi in moduli)
    print("M = {} = {}".format(prod_str, M))

    # Sous-produits Mi
    sep("Sous-produits")
    Mi_list = []
    for i in range(k):
        Mi = M // moduli[i]
        Mi_list.append(Mi)
        print("M{} = M / m{} = {}".format(i+1, i+1, Mi))

    # Inverses yi
    sep("Recherche des inverses (Mi * yi ≡ 1 [mi])")
    yi_list = []
    for i in range(k):
        Mi = Mi_list[i]
        mi = moduli[i]
        r = Mi % mi
        ok, yi = inv_mod(Mi, mi, show=False)  # on ne spamme pas Euclide ici
        if not ok:
            print("Impossible de trouver l'inverse de {} modulo {} (devrait être possible ici).".format(Mi, mi))
            return
        yi_list.append(yi)
        print("{}*y{} ≡ 1 [{}] -> {} ≡ {} [{}] -> y{} ≡ {} [{}] => y{} = {}".format(
            Mi, i+1, mi, Mi, r, mi, i+1, yi, mi, i+1, yi
        ))

    # Construction (formule CRT)
    sep("Construction (formule CRT)")
    terms_str = "  +  ".join("{} * {} * {}".format(residues[i], Mi_list[i], yi_list[i]) for i in range(k))
    print("x ≡ {}  [ {} ]".format(terms_str, M))

    # Calcul des termes
    sep("Calcul des termes")
    terms = []
    for i in range(k):
        t = residues[i] * Mi_list[i] * yi_list[i]
        terms.append(t)
        print("Terme #{} = {}*{}*{} = {}".format(i+1, residues[i], Mi_list[i], yi_list[i], t))

    # Somme et réduction
    sep("Somme")
    S = sum(terms)
    print("S = {}".format(" + ".join(str(t) for t in terms)), end="")
    print(" = {}".format(S))
    sep("Réduction")
    x0 = S % M
    print("x0 = {} % {} = {}".format(S, M, x0))

    # Solution et vérifications
    sep("Solution canonique")
    print("x ≡ {}  [ {} ]".format(x0, M))

    sep("Vérifications")
    for i in range(k):
        mi = moduli[i]
        ai = residues[i]
        print("{} % {} = {}  (attendu {}){}".format(
            x0, mi, x0 % mi, ai, "  (ok)" if (x0 % mi) == ai else "  (!!)"
        ))

    # Forme générale
    sep("Forme générale")
    print("x = {} + {}*k,  k entier".format(x0, M))

# ---------- Puissance mod m - méthode "décomposition binaire" ----------
def pow_mod_verbose(a, e, m):
    sep('Puissance mod m - méthode "décomposition binaire"')

    # Garde-fous
    if m <= 0:
        print("Le module doit être > 0")
        return
    if e < 0:
        print("Exposant négatif non pris en charge.")
        return
    if e == 0:
        print("Objectif : calculer {}^0  [{}]".format(a, m))
        print("{}^0 ≡ 1  [{}]".format(a, m))
        sep("Résultat")
        print("{}^{} (mod {}) = {}".format(a, e, m, 1 % m))
        return

    print("a = {}".format(a))
    print("e = {}".format(e))
    print("m = {}".format(m))
    print("\nObjectif : calculer {}^{}  [{}]".format(a, e, m))

    # Helpers d'affichage compacts
    def small_rep(x, mod):
        y = x % mod
        if y > mod // 2:
            return str(y - mod)  # ex: 11 -> -3 (si m=14)
        return str(y)

    def power2_repr(a_sym, k):
        if k == 0:
            return "{}".format(a_sym)
        s = "({}^2)".format(a_sym)
        for _ in range(1, k):
            s = "({}^2)".format(s)
        return s

    # 1) Écriture de l’exposant en base 2
    bits = []
    k = 0
    t = e
    while t > 0:
        if (t & 1) == 1:
            bits.append(k)
        t >>= 1
        k += 1
    bits_desc = sorted(bits, reverse=True)

    # Sommes numériques et en puissances de 2
    somme_num = " + ".join(str(1 << k) for k in bits_desc)
    somme_pow2 = " + ".join("2^{}".format(k) for k in bits_desc)
    print("\n1) Écriture de l’exposant en base 2")
    print("{} = {} = {}".format(e, somme_num, somme_pow2))

    # 2) Décomposition de la puissance
    print("\n2) Décomposition de la puissance")
    droite_pow2 = " + ".join("2^{}".format(k) for k in bits_desc)
    print("{}^{} = {}^({})".format(a, e, a, droite_pow2))
    print("     = " + " * ".join("{}^(2^{})".format(a, k) for k in bits_desc))
    print("     = " + " * ".join(power2_repr(a, k) for k in bits_desc))

    # 3) Calculs modulo m (paliers de carrés)
    print("\n3) Calculs modulo {} (paliers)".format(m))
    pow_values = {}  # k -> a^(2^k) mod m
    # k=0
    val = a % m
    pow_values[0] = val
    print("- {} ≡ {}  [ {} ]".format(a, small_rep(a, m), m))
    # k>=1
    max_k = bits_desc[0] if bits_desc else 0
    prev_val = val
    for kk in range(1, max_k + 1):
        raw_sq = (prev_val * prev_val)
        new_val = raw_sq % m
        # Affichage du palier kk en s'appuyant sur le palier kk-1 réduit
        print("- ({})^2  =>  ({}^2) = {}  =>  ≡ {}  [ {} ]".format(
            power2_repr(a, kk-1),
            small_rep(prev_val, m),
            raw_sq,
            small_rep(new_val, m),
            m
        ))
        pow_values[kk] = new_val
        prev_val = new_val

    # 4) Assemblage des facteurs utiles (bits à 1)
    print("\n4) Assemblage des facteurs utiles ({})".format(", ".join("2^{}".format(k) for k in bits_desc)))
    # Ligne symbolique
    print("{}^{} ≡ {}  [ {} ]".format(
        a, e, " * ".join(power2_repr(a, k) for k in bits_desc), m
    ))
    # Ligne numérique (représentants compacts)
    factors_str = " * ".join(small_rep(pow_values[k], m) for k in bits_desc)
    print("     ≡ {}  [ {} ]".format(factors_str, m))

    # Multiplication pas à pas
    acc = 1 % m
    if bits_desc:
        # Premier facteur
        acc = pow_values[bits_desc[0]] % m
        print("     -> {} (premier facteur)".format(small_rep(acc, m)))
        # Puis enchaînement
        for kk in bits_desc[1:]:
            before = acc
            acc = (acc * pow_values[kk]) % m
            print("     -> ({} * {}) % {} = {}".format(
                small_rep(before, m),
                small_rep(pow_values[kk], m),
                m,
                small_rep(acc, m)
            ))

    # 5) Résultat
    sep("Résultat")
    print("{}^{} (mod {}) = {}".format(a, e, m, acc))
    print("(Vérif rapide : {} % {} = {})".format(acc, m, acc % m))

def show_factorization_and_option_gcd():
    sep("Décomp. facteurs premiers")
    k = int(input("Nombre d'entiers (1 ou 2) = "))
    if k not in (1, 2):
        print("Choix invalide (1 ou 2).")
        return

    # Premier entier
    n1 = int(input("n1 = "))
    sep("n1 : échelle")
    f1 = prime_factors_ladder(n1)
    # Affichage n1 = ...
    absn1 = abs(n1)
    print("{} = {}".format(n1, _format_factorization({p:e for p,e in f1.items() if p != -1} if -1 in f1 else f1)))

    if k == 1:
        return

    # Second entier
    n2 = int(input("n2 = "))
    sep("n2 : échelle")
    f2 = prime_factors_ladder(n2)
    print("{} = {}".format(n2, _format_factorization({p:e for p,e in f2.items() if p != -1} if -1 in f2 else f2)))

    # PGCD par facteurs (exposants minimum)
    sep("PGCD par facteurs")
    # on ignore le -1 éventuel pour le PGCD
    f1pos = {p:e for p,e in f1.items() if p > 1}
    f2pos = {p:e for p,e in f2.items() if p > 1}
    common = {}
    for p in f1pos:
        if p in f2pos:
            common[p] = min(f1pos[p], f2pos[p])

    # valeur numérique du pgcd
    pgcd_val = 1
    for p, e in common.items():
        # puissance p**e
        v = 1
        for _ in range(e):
            v *= p
        pgcd_val *= v

    # Affichage type "PGCD(336,126) = 2 * 3 * 7 = 42"
    if common:
        fact_str = _format_factorization(common)
        print("PGCD({}, {}) = {} = {}".format(n1, n2, fact_str, pgcd_val))
    else:
        print("PGCD({}, {}) = 1".format(n1, n2))

# ---------- Menu "one-shot" ----------
def menu():
    sep("MENU")
    print("1) Décomp. facteurs premiers")
    print("2) PGCD / Bézout")
    print("3) Inverse mod m")
    print("4) Équations modulaires")
    print("5) Tables (Z / Z_n)")
    print("6) CRT (théorème des restes chinois)")
    print("7) Puissance mod m")
    print("8) Quitter")
    choice = input("> Choix : ").strip()

    if choice == "1":
        try:
            show_factorization_and_option_gcd()
        except:
            print("Entrée invalide.")
    elif choice == "2":
        try:
            a = int(input("a = "))
            b = int(input("b = "))
            egcd_verbose(a, b, show=True, show_back=True)
        except:
            print("Entrée invalide.")
    elif choice == "3":
        try:
            a = int(input("a = "))
            m = int(input("m = "))
            inv_mod(a, m, show=True)
        except:
            print("Entrée invalide.")
    elif choice == "4":
        equations_menu_option()
    elif choice == "5":
        try:
            sep("Tables (Z / Z_n)")
            print("Espace ?")
            print("  1 = Z (entiers)")
            print("  2 = Z_n (modulo n)")
            space = input("> Choix espace : ").strip()

            print("Opération ?")
            print("  1 = addition")
            print("  2 = multiplication")
            print("  3 = les deux")
            op_ch = input("> Choix opération : ").strip()

            def do_ops(do_add, do_mul, in_Z, in_Zn):
                if in_Z:
                    print("Intervalle en Z :")
                    s = int(input("  début = "))
                    e = int(input("  fin   = "))
                    if do_add: table_Z(s, e, "+")
                    if do_mul: table_Z(s, e, "*")
                else:
                    n = int(input("Module n (>0) : "))
                    if do_add: table_Zn(n, "+")
                    if do_mul: table_Zn(n, "*")

            do_add = (op_ch == "1" or op_ch == "3")
            do_mul = (op_ch == "2" or op_ch == "3")

            if space == "1":
                do_ops(do_add, do_mul, True, False)
            elif space == "2":
                do_ops(do_add, do_mul, False, True)
            else:
                print("Choix d'espace invalide.")
        except:
            print("Entrée invalide.")
    elif choice == "6":
        try:
            solve_system_crt_coprime()
        except:
            print("Entrée invalide.")
    elif choice == "7":
        try:
            a = int(input("a = "))
            e = int(input("e (exposant) = "))
            m = int(input("m (module > 0) = "))
            pow_mod_verbose(a, e, m)
        except:
            print("Entrée invalide.")
    elif choice == "8":
        print("Quitter le programme.")
        return
    else:
        print("Choix inconnu.")

# ---------- Lancer ----------
menu()
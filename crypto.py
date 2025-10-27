# ================================================
#  Bézout, inverse modulaire et congruences (NumWorks)
#  + Tables d'addition/multiplication en Z et Z_n (+ Z_n et Z_n*)
#  + Équation affine ax + c = 0 dans Z_n
#  + Cours (formules utiles)
#  - Euclide étendu pas à pas (traçage) + remontée
#  - Inverse mod m (avec étapes)
#  - Résolution de ax ≡ b [m] (tous cas)
#  - Système modulaire / CRT / Puissance mod m
#  Auteur : toi ;)
#  Version menu "one-shot" (pas de boucle)
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

# ---------- Équation affine ax + c = 0 dans Z_n ----------
def solve_affine_eq(a, c, n, show=True):
    sep("Équation {}x + {} = 0  [ {} ]".format(a, c, n))
    if n <= 0:
        print("Le module doit être > 0")
        return False, None, None, None
    b = (-c) % n
    print("Réécriture : {}x ≡ -{} ≡ {}  [ {} ]".format(a, c, b, n))
    return solve_congruence(a, b, n, show=show, list_rep=True)

# ---------- Résolution d'équations : wrappers pour formes générales ----------

def solve_ax_plus_c_eq_b(a, c, b, m):
    """Résout a x + c ≡ b [m] en affichant la mise en forme puis solve_congruence."""
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
    """Résout x + c ≡ b [m] (affiche la réduction élémentaire puis résout 1*x ≡ b-c [m])."""
    sep("0) Mise en forme / réduction modulo {}".format(m))
    print("{} ≡ {}  [ {} ]".format(c, c % m, m))
    print("{} ≡ {}  [ {} ]".format(b, b % m, m))
    print("=> x ≡ {}  [ {} ]".format((b - c) % m, m))
    return solve_congruence(1, (b - c) % m, m, show=True, list_rep=True)

def solve_batch_affine_zero_same_mod():
    """Résout plusieurs équations du type a_i x + c_i ≡ 0 dans Z/mZ (ex: cas de Z/31Z)."""
    sep("Résolutions dans Z/mZ")
    m = int(input("  m (module > 0) = "))
    if m <= 0:
        print("Le module doit être > 0")
        return
    k = int(input("  nombre d'équations k = "))
    if k <= 0:
        print("k doit être >= 1")
        return
    for i in range(1, k+1):
        print("Équation #{} : a x + c = 0".format(i))
        a = int(input("  a = "))
        c = int(input("  c = "))
        # On réutilise la mise en forme optimale pour ce cas
        sep("Équation {}x + {} = 0  [ {} ]".format(a, c, m))
        print("Réécriture : {}x ≡ -{} ≡ {}  [ {} ]".format(a, c, (-c) % m, m))
        solve_congruence(a, (-c) % m, m, show=True, list_rep=True)

def equations_menu_option3():
    """Sous-menu Option 3 : couvre ax≡b, ax+c≡0, ax+c≡b, x+c≡b, et série dans Z/mZ."""
    sep("Équations - choisir une forme")
    print("1) a x ≡ b  [ m ]")
    print("2) a x + c ≡ 0  [ m ]")
    print("3) a x + c ≡ b  [ m ]")
    print("4) x + c ≡ b  [ m ]")
    print("5) Plusieurs (a_i x + c_i = 0) dans Z/mZ")
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
            # Cette forme est déjà gérée par solve_affine_eq avec la bonne mise en forme
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

        elif ch == "5":
            solve_batch_affine_zero_same_mod()

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

# ---------- Fusion de deux congruences ----------
def combine_two_congruences(a1, m1, a2, m2):
    """
    Combine x ≡ a1 [m1] et x ≡ a2 [m2].
    Retourne (ok, a, M) avec x ≡ a [M], ou (False, None, None) si incompatible.
    Affiche les étapes (Euclide sur m1,m2 puis résolution).
    """
    sep("Combinaison de deux congruences")
    print("x ≡ {}  [ {} ]".format(a1, m1))
    print("x ≡ {}  [ {} ]".format(a2, m2))

    # Étape 1 : d = gcd(m1,m2)
    d, _, _ = egcd_verbose(m1, m2, show=True, show_back=False)

    # Compatibilité
    diff = (a2 - a1)
    if diff % d != 0:
        print("Incompatible : {} ne divise pas ({}-{}).".format(d, a2, a1))
        return False, None, None

    # Étape 2 : résoudre m1' * t ≡ (a2-a1)' [m2']
    m1p = m1 // d
    m2p = m2 // d
    rhs = diff // d
    sep("Réduction pour trouver t")
    print("On résout : {} * t ≡ {}  [ {} ]".format(m1p, rhs, m2p))
    ok, t0, mod_t, d_t = solve_congruence(m1p, rhs, m2p, show=True, list_rep=False)
    if not ok:
        print("Échec inattendu sur la congruence réduite.")
        return False, None, None

    # Étape 3 : solution combinée
    M = (m1 // d) * m2
    a = (a1 + m1 * t0) % M
    sep("Solution combinée")
    print("x ≡ {}  [ {} ]".format(a, M))
    return True, a, M

# ---------- (6) Système modulaire : fusion 2 à 2 ----------
def solve_system_modular():
    """
    Demande k cong. x ≡ ai [mi] et les fusionne 2 à 2 (cas général).
    """
    sep("Système modulaire (général)")
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

    # Fusion progressive
    cur_a = residues[0]
    cur_m = moduli[0]
    for i in range(1, k):
        ok, new_a, new_m = combine_two_congruences(cur_a, cur_m, residues[i], moduli[i])
        if not ok:
            print("=> Système SANS solution.")
            return
        cur_a, cur_m = new_a, new_m

    sep("Solution du système")
    print("x ≡ {}  [ {} ]".format(cur_a, cur_m))

# ---------- (7) CRT (méthode du prof — formule directe) ----------
def solve_system_crt_coprime():
    """
    Système x ≡ ai [mi] avec mi deux à deux copremiers.
    Affiche la méthode "formule directe" du prof :
      - Système
      - Produit total M
      - Sous-produits Mi = M/mi
      - Inverses yi : Mi*yi ≡ 1 [mi] (avec réduction Mi % mi)
      - Construction x ≡ Σ ai*Mi*yi [M]
      - Calcul des termes, somme S, réduction x0 = S % M
      - Vérifications et forme générale
    """
    sep("CRT (méthode du prof — formule directe)")
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

# ---------- (8) Puissance mod m — méthode "décomposition binaire" ----------
def pow_mod_verbose(a, e, m):
    """
    Calcule a^e (mod m) en affichant TOUTES les étapes :
    - écriture binaire de e
    - décomposition a^e en produit des a^(2^k) utiles
    - paliers de carrés modulo m
    - assemblage des facteurs sélectionnés (bits à 1)
    """
    sep('Puissance mod m — méthode "décomposition binaire"')

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
        """Représentant 'court' (utilise x-mod si plus petit en valeur absolue)."""
        y = x % mod
        if y > mod // 2:
            return str(y - mod)  # ex: 11 -> -3 (si m=14)
        return str(y)

    def power2_repr(a_sym, k):
        """Texte pour a^(2^k) via emboîtement (a^2)^2..."""
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
        print("- ({})^2  ⇒  ({}^2) = {}  ⇒  ≡ {}  [ {} ]".format(
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
        print("     → {} (premier facteur)".format(small_rep(acc, m)))
        # Puis enchaînement
        for kk in bits_desc[1:]:
            before = acc
            acc = (acc * pow_values[kk]) % m
            print("     → ({} * {}) % {} = {}".format(
                small_rep(before, m),
                small_rep(pow_values[kk], m),
                m,
                small_rep(acc, m)
            ))

    # 5) Résultat
    sep("Résultat")
    print("{}^{} (mod {}) = {}".format(a, e, m, acc))
    print("(Vérif rapide : {} % {} = {})".format(acc, m, acc % m))

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
    """Affiche l'échelle de divisions et renvoie le dict {p: e}."""
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

def show_factorization_and_option_gcd():
    """Menu mini: 1 ou 2 entiers; affiche l'échelle + factorisation; si 2 -> PGCD via facteurs."""
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

# ---------- Cours (formules utiles) : version NumWorks-friendly ----------
def show_course():
    sep("Cours - Choisir une fiche")
    print("1) Euclide & Bezout")
    print("2) Inverse mod m")
    print("3) ax mod m = b")
    print("4) Zn et Zn*")
    print("5) Puissances mod")
    print("6) CRT (restes chinois)")
    ch = input("> Choix : ").strip()

    sep("Fiche")
    if ch == "1":
        print("Euclide & Bezout")
        print("- gcd(a,b)=g")
        print("- il existe x,y : a*x+b*y=g")
        print("- si g=1 => a,b copremiers")
    elif ch == "2":
        print("Inverse modulo m")
        print("- existe ssi gcd(a,m)=1")
        print("- calcul: Euclide etendu")
        print("- si p premier:")
        print("  a^(p-1) mod p = 1")
        print("  a^(-1) mod p = a^(p-2)")
    elif ch == "3":
        print("Resoudre ax mod m = b")
        print("- d=gcd(a,m)")
        print("- si d ne divise pas b -> 0 sol")
        print("- sinon: a'=a/d, b'=b/d, m'=m/d")
        print("- x0 = inv(a',m')*b' mod m'")
        print("- solutions: x = x0 mod m'")
        print("- nb de classes mod m: d")
    elif ch == "4":
        print("Zn et Zn*")
        print("- Zn = {0..n-1}")
        print("- Zn* = {a in Zn : gcd(a,n)=1}")
        print("- |Zn*| = phi(n)")
        print("- si p premier: |Zp*|=p-1")
        print("- Zp est un corps")
    elif ch == "5":
        print("Puissances modulo n")
        print("- si gcd(a,n)=1 :")
        print("  a^phi(n) mod n = 1")
        print("- si p premier :")
        print("  a^(p-1) mod p = 1 (a!=0 mod p)")
    elif ch == "6":
        print("CRT (restes chinois)")
        print("- si m1,m2 copremiers :")
        print("  x=a1 (mod m1) et x=a2 (mod m2)")
        print("  => solution unique mod m1*m2")
    else:
        print("Choix inconnu.")

# ---------- Menu "one-shot" ----------
def menu():
    sep("MENU")
    print("1) Bézout / pgcd (avec étapes + remontée)")
    print("2) Inverse mod m (avec remontée)")
    print("3) Équations modulaires (sous-menu)")
    print("4) Tables (Z / Z_n)")
    print("5) Équation ax + c = 0 (Z_n)")
    print("6) Système modulaire (x ≡ a_i [m_i])")
    print("7) CRT (méthode du prof — formule directe)")
    print("8) Puissance mod m")
    print("9) Cours (formules utiles)")
    print("10) Décomp. facteurs premiers")
    print("11) Quitter")
    choice = input("> Choix : ").strip()

    if choice == "1":
        try:
            a = int(input("a = "))
            b = int(input("b = "))
            egcd_verbose(a, b, show=True, show_back=True)
        except:
            print("Entrée invalide.")
    elif choice == "2":
        try:
            a = int(input("a = "))
            m = int(input("m = "))
            inv_mod(a, m, show=True)
        except:
            print("Entrée invalide.")
    elif choice == "3":
        equations_menu_option3()
    elif choice == "4":
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
    elif choice == "5":
        try:
            print("Résoudre a x + c = 0 dans Z_n :")
            n = int(input("  n (module > 0) = "))
            a = int(input("  a = "))
            c = int(input("  c = "))
            solve_affine_eq(a, c, n, show=True)
        except:
            print("Entrée invalide.")
    elif choice == "6":
        try:
            solve_system_modular()
        except:
            print("Entrée invalide.")
    elif choice == "7":
        try:
            solve_system_crt_coprime()
        except:
            print("Entrée invalide.")
    elif choice == "8":
        try:
            a = int(input("a = "))
            e = int(input("e (exposant) = "))
            m = int(input("m (module > 0) = "))
            pow_mod_verbose(a, e, m)
        except:
            print("Entrée invalide.")
    elif choice == "9":
        show_course()
    elif choice == "10":
        try:
            show_factorization_and_option_gcd()
        except:
            print("Entrée invalide.")
    elif choice == "11":
        print("Quitter le programme.")
        return
    else:
        print("Choix inconnu.")

# ---------- Lancer ----------
menu()
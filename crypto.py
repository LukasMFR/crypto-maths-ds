# ================================================
#  Bézout, inverse modulaire et congruences (NumWorks)
#  + Tables d'addition/multiplication en Z et Z_n
#  - Euclide étendu pas à pas (traçage)
#  - Remontée (substitutions inverses) pour Bézout
#  - Inverse mod m (avec étapes)
#  - Résolution de ax ≡ b [m] (tous cas)
#  - Tables (Z / Z_n) addition & multiplication
#  Auteur : toi ;)
#  Version menu "one-shot" (pas de boucle)
# ================================================

# ---------- Outils d'affichage ----------
def sep(title=None):
    line = "-" * 38
    if title:
        print("\n" + line)
        print(title)
        print(line)
    else:
        print("\n" + line)

def show_divisions(divs):
    # divs : liste de tuples (A, B, q, r) pour "A = q*B + r"
    for (A, B, q, r) in divs:
        print("{} = {}*{} + {}".format(A, q, B, r))

def _expr_to_string(expr, R):
    # expr : dict {index -> coeff} sur les rémanents R[i]
    # R[0]=a, R[1]=b, R[2]=..., R[k]=g, R[k+1]=0
    terms = []
    for i in range(len(R)):
        c = expr.get(i, 0)
        if c != 0:
            terms.append("{}*{}".format(c, R[i]))
    if not terms:
        return "0"
    return " + ".join(terms)

def _col_width(values):
    w = 1
    for v in values:
        l = len(str(v))
        if l > w:
            w = l
    return w + 1  # petite marge

def _print_table(header_vals, row_vals, cell_fn, title):
    """Affiche un tableau avec en-tête/étiquettes de ligne.
       header_vals: liste des colonnes
       row_vals:    liste des lignes
       cell_fn(i,j): valeur à afficher
    """
    sep(title)
    # calcul largeur colonne à partir de TOUT ce qui s'affiche
    candidates = list(header_vals) + list(row_vals)
    for i in row_vals:
        for j in header_vals:
            candidates.append(cell_fn(i, j))
    w = _col_width(candidates)

    # en-tête
    first_cell = " " * w + "|"
    print(first_cell, end="")
    for j in header_vals:
        s = str(j)
        print(s.rjust(w), end="")
    print()

    # séparation
    print("-" * (w + 1 + w * len(header_vals)))

    # lignes
    for i in row_vals:
        print(str(i).rjust(w) + "|", end="")
        for j in header_vals:
            v = cell_fn(i, j)
            print(str(v).rjust(w), end="")
        print()

# ---------- gcd silencieux (pour Z_n*) ----------
def gcd(a, b):
    a = abs(a); b = abs(b)
    while b:
        a, b = b, a % b
    return a

# ---------- Euclide étendu avec traçage + remontée ----------
def egcd_verbose(a, b, show=True, show_back=True):
    """Retourne (g,x,y) avec g=gcd(a,b) et ax+by=g.
       Affiche divisions d'Euclide et, si show_back, la remontée (substitutions)."""

    A0, B0 = a, b
    divs = []      # tuples (old_r, r, q, rem) pour "old_r = q*r + rem"
    R = [a, b]     # suite des rémanents R[0]=a, R[1]=b, R[2]=rem1, ..., R[-1]=0
    Q = [None, None]  # Q[i] tel que R[i] = R[i-2] - Q[i]*R[i-1] (pour i>=2)

    # Variables pour EGCD (coefficients s,t)
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
    k = len(R) - 2  # indice du dernier rémanent non nul

    if show:
        sep("Algorithme d'Euclide ({} , {})".format(A0, B0))
        show_divisions(divs)
        print("pgcd({}, {}) = {}".format(A0, B0, g))

    # ----- Remontée : substitutions inverses pour exprimer g en fct de a,b -----
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
    """Inverse de a modulo m si gcd(a,m)=1.
       Retourne (has_inverse, inverse_mod_m)"""
    if m <= 0:
        if show:
            print("Le module doit être > 0")
        return False, None

    # On affiche Euclide + remontée pour bien montrer les étapes
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
    """Résout a x ≡ b (mod m) en affichant les étapes.
       Retourne (has_solution, x0, modulus_reduced, d)."""
    if m <= 0:
        if show:
            print("Le module doit être > 0")
        return False, None, None, None

    sep("Résolution de {} x ≡ {}  [ {} ]".format(a, b, m))
    # Étape 1 : gcd(a,m)
    d, xg, yg = egcd_verbose(a, m, show=True, show_back=False)
    if b % d != 0:
        print("Comme {} ne divise pas {}, aucune solution.".format(d, b))
        return False, None, None, d

    # Réduction
    a1 = a // d
    b1 = b // d
    m1 = m // d
    print("On réduit : a'={}, b'={}, m'={} (d = {})".format(a1, b1, m1, d))
    print("Nouvelle équation : {} x ≡ {}  [ {} ]".format(a1, b1, m1))

    # Étape 2 : inverse de a' modulo m' (avec remontée)
    ok, inv = inv_mod(a1, m1, show=True)
    if not ok:
        print("Problème inattendu : a' et m' ne sont pas copremiers.")
        return False, None, None, d

    x0 = (inv * b1) % m1
    print("Solution de base : x0 ≡ {} * {} ≡ {}  [ {} ]".format(inv, b1, x0, m1))
    print("Vérif : ({}*{}) % {} = {}  (doit ≡ {})"
          .format(a, x0, m, (a*x0) % m, b % m))

    print("Forme générale des solutions : x ≡ {}  [ {} ]".format(x0, m1))
    if d > 1:
        print("Donc modulo {}, on obtient {} solutions distinctes :".format(m, d))
        reps = []
        for k in range(d):
            reps.append((x0 + k*m1) % m)
        reps = sorted(list(set(reps)))
        if list_rep:
            print("Représentants (mod {}): {}".format(m, reps))

    return True, x0, m1, d

# ---------- Tables Z / Z_n ----------
def table_Z(start, end, op):
    """Tables en Z (sur [start..end]) pour op='+' ou '*'."""
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
    """Tables en Z_n (0..n-1) pour op='+' ou '*' + affichage Z_n et Z_n*."""
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

    # --- Ensembles Z_n et Z_n* (avec list comprehensions pour NumWorks) ---
    Zn_list = [i for i in range(n)]
    print("\nZ_{} = {{ {} }}".format(n, ", ".join([str(x) for x in Zn_list])))

    Zn_star = [a for a in range(n) if gcd(a, n) == 1]
    print("Z_{}* = {{ {} }}".format(n, ", ".join([str(x) for x in Zn_star])))

# ---------- Menu "one-shot" (pas de boucle) ----------
def menu():
    sep("MENU")
    print("1) Bézout / pgcd (avec étapes + remontée)")
    print("2) Inverse mod m (avec remontée)")
    print("3) Résoudre a x ≡ b [m]")
    print("4) Tables (Z / Z_n)")
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
        try:
            a = int(input("a = "))
            b = int(input("b = "))
            m = int(input("m = "))
            solve_congruence(a, b, m, show=True)
        except:
            print("Entrée invalide.")
    elif choice == "4":
        try:
            sep("Tables (Z / Z_n)")
            # --- PROMPTS MULTILIGNES (lisibles sur NumWorks) ---
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
    else:
        print("Choix inconnu.")

# ---------- Lancer (exécution unique) ----------
menu()
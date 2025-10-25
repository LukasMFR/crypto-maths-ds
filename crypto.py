# ================================================
#  Bézout, inverse modulaire et congruences (NumWorks)
#  + Tables d'addition/multiplication en Z et Z_n (+ Z_n et Z_n*)
#  + Équation affine ax + c = 0 dans Z_n
#  + Cours (formules utiles)
#  - Euclide étendu pas à pas (traçage) + remontée
#  - Inverse mod m (avec étapes)
#  - Résolution de ax ≡ b [m] (tous cas)
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
    print("Solution de base : x0 ≡ {} * {} ≡ {}  [ {} ]".format(inv, b1, x0, m1))
    print("Vérif : ({}*{}) % {} = {}  (doit ≡ {})".format(a, x0, m, (a*x0) % m, b % m))
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

# ---------- Équation affine ax + c = 0 dans Z_n ----------
def solve_affine_eq(a, c, n, show=True):
    sep("Équation {}x + {} = 0  [ {} ]".format(a, c, n))
    if n <= 0:
        print("Le module doit être > 0")
        return False, None, None, None
    b = (-c) % n
    print("Réécriture : {}x ≡ -{} ≡ {}  [ {} ]".format(a, c, b, n))
    return solve_congruence(a, b, n, show=show, list_rep=True)

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
    print("3) Résoudre a x ≡ b [m]")
    print("4) Tables (Z / Z_n)")
    print("5) Équation ax + c = 0 (Z_n)")
    print("6) Cours (formules utiles)")
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
        show_course()
    else:
        print("Choix inconnu.")

# ---------- Lancer ----------
menu()
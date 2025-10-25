# ================================================
#  Bézout, inverse modulaire et congruences (NumWorks)
#  - Euclide étendu pas à pas (traçage)
#  - Remontée (substitutions inverses) pour Bézout
#  - Inverse mod m (avec étapes)
#  - Résolution de ax ≡ b [m] (tous cas)
#  Auteur : toi ;)
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
    # Indice du dernier rémanent non nul (g)
    k = len(R) - 2  # car R[-1] == 0

    if show:
        sep("Algorithme d'Euclide ({} , {})".format(A0, B0))
        show_divisions(divs)
        print("pgcd({}, {}) = {}".format(A0, B0, g))

    # ----- Remontée : substitutions inverses pour exprimer g en fct de a,b -----
    if show and show_back and k >= 2:
        sep("Remontée (combinaison linéaire)")
        # On part de : R[k] = R[k-2] - Q[k]*R[k-1]
        print("{} = {} - {}*{}".format(R[k], R[k-2], Q[k], R[k-1]))
        # expr représente R[k] comme combinaison linéaire des R[i]
        expr = {k-2: 1, k-1: -Q[k]}

        # Puis on remplace R[k-1], R[k-2], ..., jusqu’à n’avoir que a=R[0] et b=R[1]
        for j in range(k-1, 1, -1):
            cj = expr.get(j, 0)
            if cj == 0:
                continue
            # R[j] = R[j-2] - Q[j]*R[j-1]
            print("Remplacer {} par {} - {}*{}".format(R[j], R[j-2], Q[j], R[j-1]))

            # Mise à jour de l'expression
            expr.pop(j, None)
            expr[j-2] = expr.get(j-2, 0) + cj
            expr[j-1] = expr.get(j-1, 0) - cj * Q[j]

            # Afficher l'état courant
            print("=> {} = {}".format(R[k], _expr_to_string(expr, R)))

        # À la fin, expr = {0:x, 1:y}
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

    # Forme générale
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

# ---------- Menu interactif (sans démo) ----------
def menu():
    while True:
        sep("MENU")
        print("1) Bézout / pgcd (avec étapes + remontée)")
        print("2) Inverse mod m (avec remontée)")
        print("3) Résoudre a x ≡ b [m]")
        print("0) Quitter")
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
        elif choice == "0":
            print("Bye.")
            break
        else:
            print("Choix inconnu.")

# ---------- Lancer ----------
menu()
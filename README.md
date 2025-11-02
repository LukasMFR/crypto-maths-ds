# Bézout, inverse modulaire et congruences (NumWorks)

Programme Python complet et optimisé pour **NumWorks**, dédié à l’étude des **équations et systèmes modulaires** en arithmétique.  
Ce script permet de manipuler les congruences, le PGCD, les inverses modulaires, les puissances mod m, ainsi que le théorème des restes chinois (CRT).

> 🧮 Version : **menu “one-shot”** — sans boucle continue, optimisée pour la mémoire RAM de la NumWorks (~32 Ko max).

---

## 📚 Sommaire
1. [Présentation](#-présentation)
2. [Fonctionnalités](#-fonctionnalités)
3. [Structure du menu](#-structure-du-menu)
4. [Méthodes utilisées](#-méthodes-utilisées)
5. [Compatibilité NumWorks](#-compatibilité-numworks)
6. [Exemples d’utilisation](#-exemples-dutilisation)
7. [Auteur](#-auteur)

---

## 🧩 Présentation

Ce programme regroupe dans un même fichier Python plusieurs outils de base en **arithmétique modulaire** :

- Calcul du **PGCD et des coefficients de Bézout**
- Recherche d’**inverses modulaires**
- Résolution d’**équations linéaires modulaires**
- Génération de **tables en Z et Zₙ**
- Application du **théorème des restes chinois (CRT)** pour des systèmes à modules copremiers
- Calcul de **puissances modulo m**
- **Décomposition en facteurs premiers** (méthode de l’échelle)

Le script est conçu pour être **auto-suffisant**, pédagogique, et fonctionner sans dépendance externe.

---

## ⚙️ Fonctionnalités

| # | Fonction | Description | Exemple d’utilisation |
|---|-----------|--------------|------------------------|
| 1 | **Décomposition en facteurs premiers** | Factorise un entier et affiche les étapes de division. | `336 → 2⁴ × 3 × 7` |
| 2 | **PGCD / Bézout** | Applique l’algorithme d’Euclide étendu et montre la remontée pas à pas. | `pgcd(252, 198)` |
| 3 | **Inverse mod m** | Calcule et vérifie l’inverse modulaire si gcd(a, m) = 1. | `a = 7, m = 26 → 15` |
| 4 | **Équations modulaires** | Résout `a x ≡ b [m]`, `a x + c ≡ b [m]`, `x + c ≡ b [m]` avec mise en forme automatique. | `4x + 5 ≡ 0 [21]` |
| 5 | **Tables Z / Zₙ** | Affiche la table d’addition ou de multiplication sur un intervalle ou modulo n. | `Z₇* = {1,2,3,4,5,6}` |
| 6 | **CRT (formule directe)** | Résout un système de congruences copremières avec la formule du théorème des restes chinois. | `x ≡ 7 [10], x ≡ 3 [7], x ≡ 6 [13] → x ≡ 227 [910]` |
| 7 | **Puissance mod m** | Calcule `a^e mod m` étape par étape (méthode de décomposition binaire). | `3¹³ mod 17` |

---

## 🧮 Structure du menu

Lors de l’exécution, le programme affiche :

```

--- MENU ---

1. Décomp. facteurs premiers
2. PGCD / Bézout
3. Inverse mod m
4. Équations modulaires
5. Tables (Z / Z_n)
6. CRT (formule directe)
7. Puissance mod m
8. Quitter

````

Chaque option exécute une fonction spécifique.  
Le script ne relance pas le menu après une opération (version “one-shot”), ce qui permet d’économiser de la mémoire RAM sur la NumWorks.

---

## 🧠 Méthodes utilisées

- **Algorithme d’Euclide étendu** pour le PGCD et les coefficients de Bézout.
- **Résolution d’équations modulaires** via l’inverse modulaire.
- **Théorème des restes chinois (CRT)** pour les systèmes à modules copremiers :
  \[
  x ≡ a_i \pmod{m_i}, \quad \text{avec } \gcd(m_i, m_j) = 1
  \]
  et construction :
  \[
  x ≡ \sum a_i M_i y_i \pmod{M}
  \]
  où \( M_i = M/m_i \) et \( y_i = M_i^{-1} \pmod{m_i} \)
- **Méthode de décomposition binaire** pour le calcul de puissances modulaires rapides.

---

## 🧩 Compatibilité NumWorks

- ✅ Compatible **Epsilon / Upsilon Custom / Omega / Khi**
- ✅ Testé sur modèle **N0110**
- ⚙️ Taille maximale recommandée : **~30 Ko**
- ⚠️ Pas de boucle principale (`while True:`) pour préserver la mémoire.

### 📦 Import
Copiez le fichier dans le dossier Python de votre NumWorks (ou via Epsilon Online)  
et exécutez simplement :
```python
menu()
````

---

## 🔢 Exemple d’utilisation : CRT

```
> Choix : 6
--- CRT (formule directe) ---
Nombre d'équations k = 3
Equation #1 :
  a1 = 7
  m1 (>0) = 10
Equation #2 :
  a2 = 3
  m2 (>0) = 7
Equation #3 :
  a3 = 6
  m3 (>0) = 13
--- Produit total ---
M = 10 * 7 * 13 = 910
--- Inverses ---
y1 = 1, y2 = 2, y3 = 8
--- Construction ---
x ≡ 7*91*1 + 3*130*2 + 6*70*8 [910]
x0 = 227
✅ Résultat : x ≡ 227 [910]
```

---

## 👨‍💻 Auteur

**Lukas Mauffré**
Étudiant à l’ECE Paris — Passionné de mathématiques, d’algorithmique et de NumWorks.
Version : *menu “one-shot”, allégée pour NumWorks (RAM 32 Ko)*

---

> 💡 *Dernière mise à jour : novembre 2025*

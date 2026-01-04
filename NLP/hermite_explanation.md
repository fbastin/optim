# Explication Détaillée : Construction d'une Fonction C¹ avec Toutes les Propriétés Requises

## Objectif

Construire une fonction univariée continuellement différentiable (classe C¹) possédant :
1. Un **point selle** (ni minimiseur ni maximiseur)
2. Au moins un **minimiseur local** (non-global)
3. Au moins un **maximiseur local**
4. Un **minimiseur global**

## Stratégie Générale

La fonction est construite par **morceaux** en trois zones :
- **Zone 1** (x ≤ -1) : Polynômes d'Hermite cubiques
- **Zone 2** (-1 ≤ x ≤ 1) : Polynôme de degré 5 avec point selle
- **Zone 3** (x ≥ 1) : Polynômes d'Hermite cubiques

Cette approche hybride permet de :
- Garantir mathématiquement un point selle en x = 0
- Assurer la continuité C¹ aux jonctions x = -1 et x = 1
- Contrôler précisément la position des extrema

---

## 1. Polynômes d'Hermite Cubiques

### Principe

Un polynôme d'Hermite cubique entre deux points (x₀, y₀) et (x₁, y₁) est défini par :

```
H(x) = h₀₀(t)·y₀ + h₁₀(t)·dx·dy₀ + h₀₁(t)·y₁ + h₁₁(t)·dx·dy₁
```

où :
- t = (x - x₀)/(x₁ - x₀) ∈ [0,1] (paramètre normalisé)
- dx = x₁ - x₀ (largeur de l'intervalle)
- dy₀, dy₁ sont les dérivées spécifiées aux extrémités
- h₀₀, h₁₀, h₀₁, h₁₁ sont les fonctions de base d'Hermite

### Fonctions de Base

Les fonctions de base d'Hermite garantissent les conditions d'interpolation :

```
h₀₀(t) = 2t³ - 3t² + 1       → H(x₀) = y₀
h₁₀(t) = t³ - 2t² + t         → H'(x₀) = dy₀
h₀₁(t) = -2t³ + 3t²           → H(x₁) = y₁
h₁₁(t) = t³ - t²              → H'(x₁) = dy₁
```

### Propriété Clé : Continuité C¹

Entre deux segments consécutifs [xᵢ, xᵢ₊₁] et [xᵢ₊₁, xᵢ₊₂], si on impose que la dérivée à droite du premier segment égale la dérivée à gauche du second segment (dyᵢ₊₁), alors **la fonction globale est automatiquement C¹**.

### Implémentation

```julia
function hermite_cubic(x, x0, x1, y0, y1, dy0, dy1)
    t = (x - x0) / (x1 - x0)
    h00 = 2*t^3 - 3*t^2 + 1
    h10 = t^3 - 2*t^2 + t
    h01 = -2*t^3 + 3*t^2
    h11 = t^3 - t^2
    
    dx = x1 - x0
    return h00*y0 + h10*dx*dy0 + h01*y1 + h11*dx*dy1
end
```

---

## 2. Zone Centrale : Polynôme avec Point Selle

### Construction du Point Selle

Pour avoir un **point selle** en x = 0, il faut que :
- f'(0) = 0 (point critique)
- f''(0) = 0 (point d'inflexion)

On utilise un polynôme de la forme :

```
f(x) = ax⁵ + bx³ + d
```

Ses dérivées sont :
```
f'(x) = 5ax⁴ + 3bx²
f''(x) = 20ax³ + 6bx
```

En x = 0 :
```
f'(0) = 0  ✓  (automatiquement satisfait)
f''(0) = 0  ✓  (automatiquement satisfait)
```

### Raccordement C¹ aux Jonctions

Le polynôme doit se raccorder de manière C¹ avec les zones 1 et 3 aux points x = -1 et x = 1.

**Conditions à satisfaire :**
1. f(-1) = y_{zone1}(-1) (continuité de f)
2. f(1) = y_{zone3}(1) (continuité de f)
3. f'(-1) = dy_{zone1}(-1) (continuité de f')
4. f'(1) = dy_{zone3}(1) (continuité de f')

**Système d'équations :**

À partir de f(x) = ax⁵ + bx³ + d :

```
f(-1) = -a - b + d = y₁
f(1) = a + b + d = y₃
f'(-1) = 5a + 3b = dy₁
f'(1) = 5a + 3b = dy₃
```

Notez que f'(-1) = f'(1) par symétrie du polynôme. On prend donc la moyenne :

```
dy_moy = (dy₁ + dy₃)/2
```

**Résolution :**

Addition de f(-1) et f(1) :
```
2d = y₁ + y₃  →  d = (y₁ + y₃)/2
```

Soustraction :
```
2a + 2b = y₃ - y₁  →  a + b = (y₃ - y₁)/2
```

De la condition sur les dérivées :
```
5a + 3b = dy_moy
```

En combinant les deux dernières équations :
```
b = (y₃ - y₁)/2 - a
5a + 3((y₃ - y₁)/2 - a) = dy_moy
5a + 3(y₃ - y₁)/2 - 3a = dy_moy
2a = dy_moy - 3(y₃ - y₁)/2

→ a = (dy_moy - 3(y₃ - y₁)/2)/2
```

Puis :
```
b = (y₃ - y₁)/2 - a
```

### Code

```julia
a_poly = (avg_dy - 3*(y_at_plus1 - y_at_minus1)/2) / 2
b_poly = (y_at_plus1 - y_at_minus1)/2 - a_poly
d_poly = (y_at_minus1 + y_at_plus1) / 2
```

---

## 3. Contrôle des Extrema

### Minimiseur Local (Zone 1)

Dans la zone 1, on place un point avec **dy = 0** (dérivée nulle) :
```julia
x1_points = [-4.0, -3.0, -2.0, -1.0]
y1_points = [-1.0, -2.5, -3.0, -2.0]
dy1_points = [-1.5, -0.8, 0.0, 1.5]  # dy = 0 en x = -2
```

Le polynôme d'Hermite garantit que f'(-2) = 0, créant un **minimiseur local**.

### Maximiseur Local (Zone 3)

Dans la zone 3, on peut créer un maximiseur en plaçant des points qui forment une "bosse" :
```julia
x3_points = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
y3_points = [2.0, 3.0, 1.5, -1.0, -4.5, -3.0, -1.5]
```

La forme naturelle de l'interpolation d'Hermite crée un maximum local autour de x ≈ 2.

### Minimiseur Global (Zone 3)

On place un point bas avec **dy = 0** :
```julia
dy3_points = [1.5, 0.5, -2.0, -3.0, 0.0, 2.5, 1.5]  # dy = 0 en x = 5
```

Comme y(5) = -4.5 est plus bas que tous les autres minimiseurs locaux, c'est le **minimiseur global**.

---

## 4. Recherche Automatique des Points Critiques

Le code recherche numériquement les points où f'(x) ≈ 0 :

```julia
for i in 10:length(x_range)-10
    if abs(fp_vals[i]) < tolerance
        # Vérifier que c'est un minimum local de |f'|
        is_local_min = true
        for j in -5:5
            if abs(fp_vals[i]) > abs(fp_vals[i+j])
                is_local_min = false
                break
            end
        end
        
        if is_local_min
            # Classifier avec f''(x)
            fppc = f_double_prime(x_crit)
            if fppc > 0
                type_pt = "Minimiseur local"
            elseif fppc < 0
                type_pt = "Maximiseur local"
            else
                type_pt = "Point selle"
            end
        end
    end
end
```

### Classification

Pour chaque point critique (f'(x) ≈ 0), on calcule f''(x) :
- **f''(x) > 0** → Minimiseur local (concavité vers le haut)
- **f''(x) < 0** → Maximiseur local (concavité vers le bas)
- **f''(x) ≈ 0** → Point selle (point d'inflexion)

---

## 5. Vérification de la Continuité C¹

Le code vérifie explicitement aux jonctions x = -1 et x = 1 :

```julia
f_gauche = f(xj - ε)
f_droite = f(xj + ε)
fp_gauche = f_prime(xj - ε)
fp_droite = f_prime(xj + ε)

delta_f = |f_droite - f_gauche|   # Doit être ≈ 0
delta_fp = |fp_droite - fp_gauche| # Doit être ≈ 0
```

Si delta_f < 10⁻⁶ et delta_fp < 10⁻⁴, la fonction est C¹ à cette jonction.

---

## 6. Visualisation et Export

Le graphe affiche :
- La fonction f(x) en bleu
- Les points de contrôle en bleu clair
- Les jonctions x = -1 et x = 1 en lignes pointillées grises
- Le point selle en orange (◆)
- Les minimiseurs locaux en vert (●)
- Le minimiseur global en vert foncé (★)
- Les maximiseurs locaux en rouge (●)

Le graphe est exporté en PDF avec :
```julia
savefig(plt, "fonction_hermite.pdf")
```

---

## Résumé des Garanties Mathématiques

1. **Point selle garanti** : Le polynôme ax⁵ + bx³ + d satisfait f'(0) = 0 et f''(0) = 0 par construction

2. **Continuité C¹ garantie** : Les polynômes d'Hermite cubiques assurent automatiquement la continuité de f et f' entre segments consécutifs

3. **Extrema contrôlés** : En spécifiant dy = 0 aux points désirés, on force les extrema locaux

4. **Raccordement exact** : Les coefficients a, b, d sont calculés pour que le polynôme central se raccorde exactement avec les zones latérales

Cette méthode combine rigueur mathématique et flexibilité pour obtenir exactement les propriétés souhaitées.
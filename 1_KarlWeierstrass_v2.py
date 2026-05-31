"""
Moj doktorski rad LOTO 7/39 predikcije: 
- visoka nauka 
- teska matematika 
- briljantno programiranje
"""




"""
Karl Weierstrass algorithm for Loto 7/39 prediction — v2


1. Weierstrass-ova funkcija nad svih 4624 do sad izvucenih kombinacija.       
2. Aparati:  
    a) Brownovo kretanje, 
    b) Hurst eksponent (R/S analiza), 
    c) Fraktalna dimenzija.      
3. nad svakim aparatom svaki Test:  
    a) Hurst eksponent,   
    b) Autokorelacija (ACF),  
    c) Mutual information,  
    d) Sample / permutation entropy,  
    e) NIST baterija testova nasumičnosti. 

3x5 = 15 modela.            

lex-indeks = leksikografski indeks

Strukturu testiram nad nizom od 4624 lex-indeksa (iz skalar 1..15.380.937).      
Analiza rezultata i izbor sta je najbolje za predikciju sledece loto kombinacije. 

Niz nije 7 kolona nego jedan skalarni niz: svaka Loto 7/39 kombinacija se mapira u lex-indeks 1..C(39,7).
Nad tim nizom od 4624 lex-indeksa posmatram "Weierstrass/fraktalnu" krivu: 
možda nema tangente, ali možda ima skrivenu strukturu. 
Cilj nije odmah predikcija, nego prvo analiza: da li postoji prediktivna struktura.
Aparati: Brownovo kretanje, Hurst/R-S, fraktalna dimenzija.
Testovi: Hurst, ACF, mutual information, sample/permutation entropy, NIST randomness.

Niz istorijskih izvlacenja od (svih do sada izvucenih) 4624 kombinacije. 
Ceo prostor kombinacija 39C7.
Koliko od 39C7 kombinacija je istorija pogodila bar jednom (4624 / 15.380.937 ≈ 0.03%). 
Kolmogorov-Smirnov / chi-kvadrat test uniformnosti. 
Svaka izvučena kombinacija se mapira na tačan red/indeks u 39C7 prostoru.
Dobija se niz od 4624 lex-indeksa: to su stvarne 4624 tačke Karl Weierstrass krive.

Ne radi se Karl Weierstrass kriva od svih 15.380.937 tačaka 
(to nije Weierstrass kriva, to je dijagonala). 
Dobija se savršena uniformna šumna sekvencu (čist šum — nema strukture).
Nije vremenska serija — to je veličina prostora stanja. 
Fraktalnost je svojstvo putanje kroz vreme, ne prostora mogućnosti.

Radi se Weierstrass kriva od 4624 stvarno posećene tačke, 
gde je svaka tačka njen indeks u ukupnom prostoru.
Dakle f(t) je: t = redni broj izvlačenja 1..4624 f(t) = indeks te kombinacije u 39C7 prostoru.

To je osnova za Brown/Hurst/fraktal/NIST analizu.
 

KORAK 1:           Weierstrass-ova funkcija nad svih 4624 do sad izvucenih kombinacija
                   - mapiranje izvucenih kombinacija u lex-indeks iz prostora 39C7
                   - formiranje krive f(t)
                   - vizualizacija krive (PNG)
                   - sazetak (TXT)

KORAK 2a:          Aparat - Brownovo kretanje
                   - inkrementi dX(t) = f(t+1) - f(t)
                   - stacionarnost (prva polovina vs druga polovina)
                   - normalnost (Shapiro-Wilk, Anderson-Darling)
                   - varijansa scaling: Var(X(t+k)-X(t)) vs k
                     (Brown: linearno, slope ~ 1 u log-log)

KORAK 2b:          Aparat - Hurst eksponent (R/S analiza)
                   - R/S analiza nad f(t)
                   - poredjenje sa Brown/random referencom H ~ 0.5
                   - log-log grafikon R/S vs window

KORAK 2c:          Aparat - Fraktalna dimenzija
                   - Higuchi fractal dimension nad f(t)
                   - Katz fractal dimension nad f(t)
                   - poredjenje glatko/random/fraktalno

KORAK 2a3a:          Aparat 2a Brownovo kretanje + Test 3a Hurst eksponent
                     - Hurst nad Brown-putanjom iz centriranih inkremenata
                     - Hurst nad samim inkrementima
                     - shuffled Brown referenca

KORAK 2a3b:          Aparat 2a Brownovo kretanje + Test 3b Autokorelacija (ACF)
                     - ACF nad centriranim Brown inkrementima
                     - ACF nad Brown-putanjom kao kontrola
                     - Ljung-Box aproksimacija i shuffled referenca

KORAK 2a3c:          Aparat 2a Brownovo kretanje + Test 3c Mutual Information
                     - MI nad centriranim Brown inkrementima
                     - MI nad Brown-putanjom kao kontrola
                     - shuffled MI referenca

KORAK 2a3d:          Aparat 2a Brownovo kretanje + Test 3d Sample / Permutation entropy
                     - sample entropy nad centriranim Brown inkrementima
                     - permutation entropy nad centriranim Brown inkrementima
                     - shuffled entropy referenca

KORAK 2a3e:          Aparat 2a Brownovo kretanje + Test 3e NIST baterija
                     - monobit, runs, block frequency
                     - cumulative sums
                     - approximate entropy

KORAK 2b3a:          Aparat 2b Hurst/R-S + Test 3a Hurst eksponent
                     - globalni Hurst nad f(t)
                     - rolling/local Hurst kroz vreme
                     - shuffled Hurst referenca

KORAK 2b3b:          Aparat 2b Hurst/R-S + Test 3b Autokorelacija (ACF)
                     - ACF nad rolling/local Hurst nizom
                     - ACF nad f(t) kao kontrola
                     - Ljung-Box i shuffled ACF referenca

KORAK 2b3c:          Aparat 2b Hurst/R-S + Test 3c Mutual Information
                     - MI nad rolling/local Hurst nizom
                     - MI nad f(t) kao kontrola
                     - shuffled MI referenca

KORAK 2b3d:          Aparat 2b Hurst/R-S + Test 3d Sample / Permutation entropy
                     - sample entropy nad rolling/local Hurst nizom
                     - permutation entropy nad rolling/local Hurst nizom
                     - shuffled entropy referenca

KORAK 2b3e:          Aparat 2b Hurst/R-S + Test 3e NIST baterija
                     - NIST-style testovi nad rolling H > 0.5 bitovima
                     - f(t) NIST kontrola
                     - monobit, runs, block frequency, cumulative sums, approximate entropy

KORAK 2c3a:          Aparat 2c Fraktalna dimenzija + Test 3a Hurst eksponent
                     - veza Higuchi FD i Hurst reference D ~= 2 - H
                     - rolling/local FD kroz vreme
                     - shuffled FD referenca

KORAK 2c3b:          Aparat 2c Fraktalna dimenzija + Test 3b Autokorelacija (ACF)
                     - ACF nad rolling/local FD nizom
                     - ACF nad f(t) kao kontrola
                     - Ljung-Box i shuffled ACF referenca

KORAK 2c3c:          Aparat 2c Fraktalna dimenzija + Test 3c Mutual Information
                     - MI nad rolling/local FD nizom
                     - MI nad f(t) kao kontrola
                     - shuffled MI referenca

KORAK 2c3d:          Aparat 2c Fraktalna dimenzija + Test 3d Sample / Permutation entropy
                     - sample entropy nad rolling/local FD nizom
                     - permutation entropy nad rolling/local FD nizom
                     - shuffled entropy referenca

KORAK 2c3e:          Aparat 2c Fraktalna dimenzija + Test 3e NIST baterija
                     - NIST-style testovi nad rolling FD > median(FD) bitovima
                     - f(t) NIST kontrola
                    - monobit, runs, block frequency, cumulative sums, approximate entropy
"""


import csv
import math
import os
import time
from datetime import timedelta

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


T0 = time.time()

CSV_DRAWS = "/data/loto7_4624_k43.csv"
CSV_ALL_COMBOS = "/data/kombinacije_39C7.csv"

HERE = os.path.dirname(os.path.abspath(__file__))
PNG_PATH = os.path.join(HERE, "1_KarlWeierstrass_v2_1.png")
PNG_PATH_2A = os.path.join(HERE, "1_KarlWeierstrass_v2_2a.png")
PNG_PATH_2B = os.path.join(HERE, "1_KarlWeierstrass_v2_2b.png")
PNG_PATH_2C = os.path.join(HERE, "1_KarlWeierstrass_v2_2c.png")
PNG_PATH_2A3A = os.path.join(HERE, "1_KarlWeierstrass_v2_2a3a.png")
PNG_PATH_2A3B = os.path.join(HERE, "1_KarlWeierstrass_v2_2a3b.png")
PNG_PATH_2A3C = os.path.join(HERE, "1_KarlWeierstrass_v2_2a3c.png")
PNG_PATH_2A3D = os.path.join(HERE, "1_KarlWeierstrass_v2_2a3d.png")
PNG_PATH_2A3E = os.path.join(HERE, "1_KarlWeierstrass_v2_2a3e.png")
PNG_PATH_2B3A = os.path.join(HERE, "1_KarlWeierstrass_v2_2b3a.png")
PNG_PATH_2B3B = os.path.join(HERE, "1_KarlWeierstrass_v2_2b3b.png")
PNG_PATH_2B3C = os.path.join(HERE, "1_KarlWeierstrass_v2_2b3c.png")
PNG_PATH_2B3D = os.path.join(HERE, "1_KarlWeierstrass_v2_2b3d.png")
PNG_PATH_2B3E = os.path.join(HERE, "1_KarlWeierstrass_v2_2b3e.png")
PNG_PATH_2C3A = os.path.join(HERE, "1_KarlWeierstrass_v2_2c3a.png")
PNG_PATH_2C3B = os.path.join(HERE, "1_KarlWeierstrass_v2_2c3b.png")
PNG_PATH_2C3C = os.path.join(HERE, "1_KarlWeierstrass_v2_2c3c.png")
PNG_PATH_2C3D = os.path.join(HERE, "1_KarlWeierstrass_v2_2c3d.png")
PNG_PATH_2C3E = os.path.join(HERE, "1_KarlWeierstrass_v2_2c3e.png")
TXT_PATH = os.path.join(HERE, "1_KarlWeierstrass_v2.txt")

N_MAX = 39
K_PICK = 7
TOTAL_COMBOS = math.comb(N_MAX, K_PICK)


def read_loto_csv(path):
    rows = []
    with open(path, "r", encoding="utf-8", newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) < K_PICK:
                continue
            try:
                nums = tuple(sorted(int(x) for x in row[:K_PICK]))
            except ValueError:
                continue
            if len(nums) == K_PICK and len(set(nums)) == K_PICK:
                rows.append(nums)
    return rows


def lex_rank_1based(combo, n=N_MAX, k=K_PICK):
    """1-based lex indeks (poklapa se sa rednim brojem u kombinacije_39C7.csv)."""
    combo = tuple(sorted(combo))
    rank0 = 0
    prev = 0
    for i, value in enumerate(combo):
        remaining = k - i - 1
        for candidate in range(prev + 1, value):
            rank0 += math.comb(n - candidate, remaining)
        prev = value
    return rank0 + 1


def hurst_rs(series, min_window=8, max_window=None):
    """R/S Hurst procena: slope log(R/S) prema log(window)."""
    x = np.asarray(series, dtype=float)
    n = len(x)
    if max_window is None:
        max_window = max(min_window * 2, n // 4)

    windows = []
    w = min_window
    while w <= max_window:
        windows.append(w)
        w = int(w * 1.45) + 1

    used_windows = []
    rs_values = []
    for w in windows:
        chunks = n // w
        if chunks < 2:
            continue
        vals = []
        for i in range(chunks):
            seg = x[i * w:(i + 1) * w]
            y = seg - seg.mean()
            z = np.cumsum(y)
            r = z.max() - z.min()
            s = seg.std(ddof=1)
            if s > 0:
                vals.append(r / s)
        if vals:
            used_windows.append(w)
            rs_values.append(float(np.mean(vals)))

    used_windows = np.asarray(used_windows, dtype=float)
    rs_values = np.asarray(rs_values, dtype=float)
    slope, intercept = np.polyfit(np.log(used_windows), np.log(rs_values), 1)
    fit = intercept + slope * np.log(used_windows)
    ss_res = float(np.sum((np.log(rs_values) - fit) ** 2))
    ss_tot = float(np.sum((np.log(rs_values) - np.log(rs_values).mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(slope), float(intercept), float(r2), used_windows, rs_values


def normalize01(x):
    x = np.asarray(x, dtype=float)
    span = float(x.max() - x.min())
    return (x - x.min()) / (span + 1e-12)


def higuchi_fd(series, kmax=64):
    """Higuchi fractal dimension za 1D vremenski niz."""
    x = np.asarray(series, dtype=float)
    n = len(x)
    ks = np.arange(1, min(kmax, n // 2) + 1, dtype=int)
    lk = []
    used = []

    for k in ks:
        lm = []
        for m in range(k):
            idx = np.arange(m, n, k)
            if len(idx) < 2:
                continue
            dist = np.abs(np.diff(x[idx])).sum()
            norm = (n - 1) / ((len(idx) - 1) * k)
            lm.append((dist * norm) / k)
        if lm:
            used.append(k)
            lk.append(float(np.mean(lm)))

    used = np.asarray(used, dtype=float)
    lk = np.asarray(lk, dtype=float)
    slope, intercept = np.polyfit(np.log(1.0 / used), np.log(lk), 1)
    fit = intercept + slope * np.log(1.0 / used)
    ss_res = float(np.sum((np.log(lk) - fit) ** 2))
    ss_tot = float(np.sum((np.log(lk) - np.log(lk).mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(slope), float(intercept), float(r2), used, lk


def katz_fd(series):
    """Katz fractal dimension; racuna se nad normalizovanim f(t)."""
    x = normalize01(series)
    pts = np.column_stack((np.linspace(0.0, 1.0, len(x)), x))
    step_lengths = np.sqrt(np.sum(np.diff(pts, axis=0) ** 2, axis=1))
    total_len = float(step_lengths.sum())
    d = float(np.max(np.sqrt(np.sum((pts - pts[0]) ** 2, axis=1))))
    n = len(x)
    if total_len <= 0 or d <= 0:
        return float("nan")
    return float(np.log10(n) / (np.log10(d / total_len) + np.log10(n)))


def autocorr_values(series, max_lag=60):
    """ACF za lagove 0..max_lag."""
    x = np.asarray(series, dtype=float)
    x = x - x.mean()
    denom = float(np.dot(x, x))
    vals = []
    for lag in range(max_lag + 1):
        if lag == 0:
            vals.append(1.0)
        else:
            vals.append(float(np.dot(x[:-lag], x[lag:]) / denom))
    return np.asarray(vals, dtype=float)


def ljung_box_approx(acf_vals, n, h):
    """Ljung-Box Q aproksimacija za prvih h ACF lagova."""
    lags = np.arange(1, h + 1, dtype=float)
    q = n * (n + 2) * float(np.sum((acf_vals[1:h + 1] ** 2) / (n - lags)))
    p = float(stats.chi2.sf(q, h))
    return float(q), p


def quantile_symbols(series, bins=16):
    """Pretvara niz u diskretne simbole preko kvantil-binova."""
    x = np.asarray(series, dtype=float)
    edges = np.quantile(x, np.linspace(0.0, 1.0, bins + 1))
    edges = np.unique(edges)
    if len(edges) <= 2:
        edges = np.linspace(float(x.min()), float(x.max()), bins + 1)
    return np.digitize(x, edges[1:-1], right=False)


def discrete_mutual_information(sym_a, sym_b):
    """Diskretna mutual information u bitovima."""
    a = np.asarray(sym_a, dtype=int)
    b = np.asarray(sym_b, dtype=int)
    n = len(a)
    size_a = int(a.max()) + 1
    size_b = int(b.max()) + 1
    table = np.zeros((size_a, size_b), dtype=float)
    np.add.at(table, (a, b), 1.0)
    pxy = table / n
    px = pxy.sum(axis=1, keepdims=True)
    py = pxy.sum(axis=0, keepdims=True)
    expected = px @ py
    mask = pxy > 0
    return float(np.sum(pxy[mask] * np.log2(pxy[mask] / expected[mask])))


def mutual_information_lags(series, max_lag=60, bins=16):
    """MI(series[t], series[t+lag]) za lagove 1..max_lag."""
    symbols = quantile_symbols(series, bins=bins)
    vals = []
    for lag in range(1, max_lag + 1):
        vals.append(discrete_mutual_information(symbols[:-lag], symbols[lag:]))
    return np.asarray(vals, dtype=float)


def sample_entropy(series, m=2, r=None, max_points=1200):
    """Sample entropy; koristi poduzorak ako je niz dug radi brzine."""
    x = np.asarray(series, dtype=float)
    if len(x) > max_points:
        idx = np.linspace(0, len(x) - 1, max_points).astype(int)
        x = x[idx]
    x = (x - x.mean()) / (x.std() + 1e-12)
    if r is None:
        r = 0.2

    def _count(mm):
        templates = np.array([x[i:i + mm] for i in range(len(x) - mm + 1)])
        count = 0
        for i in range(len(templates) - 1):
            dist = np.max(np.abs(templates[i + 1:] - templates[i]), axis=1)
            count += int(np.sum(dist <= r))
        return count

    b = _count(m)
    a = _count(m + 1)
    if a == 0 or b == 0:
        return float("inf"), a, b, len(x)
    return float(-np.log(a / b)), a, b, len(x)


def permutation_entropy(series, order=3, delay=1):
    """Normalizovana permutation entropy u opsegu 0..1."""
    x = np.asarray(series, dtype=float)
    n_patterns = len(x) - delay * (order - 1)
    if n_patterns <= 0:
        return float("nan"), float("nan"), 0

    counts = {}
    for i in range(n_patterns):
        window = x[i:i + delay * order:delay]
        pattern = tuple(np.argsort(window, kind="mergesort"))
        counts[pattern] = counts.get(pattern, 0) + 1

    probs = np.asarray(list(counts.values()), dtype=float)
    probs = probs / probs.sum()
    pe = float(-np.sum(probs * np.log2(probs)))
    pe_norm = pe / np.log2(math.factorial(order))
    return pe, float(pe_norm), len(counts)


def nist_bits_from_series(series):
    """Binarizacija za NIST-style test: 1 ako je inkrement iznad medijane."""
    x = np.asarray(series, dtype=float)
    med = float(np.median(x))
    bits = (x > med).astype(int)
    if bits.sum() == 0 or bits.sum() == len(bits):
        bits = (x > x.mean()).astype(int)
    return bits


def nist_monobit(bits):
    n = len(bits)
    s_obs = abs(int(np.sum(2 * bits - 1))) / np.sqrt(n)
    p = math.erfc(s_obs / np.sqrt(2))
    return float(p), float(s_obs)


def nist_runs(bits):
    n = len(bits)
    pi = float(bits.mean())
    if abs(pi - 0.5) >= 2 / np.sqrt(n):
        return 0.0, 0, pi
    runs = int(1 + np.sum(bits[1:] != bits[:-1]))
    denom = 2 * np.sqrt(2 * n) * pi * (1 - pi)
    p = math.erfc(abs(runs - 2 * n * pi * (1 - pi)) / denom)
    return float(p), runs, pi


def nist_block_frequency(bits, block_size=128):
    n = len(bits)
    n_blocks = n // block_size
    if n_blocks == 0:
        return float("nan"), float("nan"), 0
    trimmed = bits[:n_blocks * block_size].reshape(n_blocks, block_size)
    props = trimmed.mean(axis=1)
    chi2 = 4 * block_size * float(np.sum((props - 0.5) ** 2))
    p = float(stats.chi2.sf(chi2, n_blocks))
    return p, chi2, n_blocks


def nist_cumulative_sums(bits):
    x = 2 * bits - 1
    walk = np.cumsum(x)
    z = int(np.max(np.abs(walk)))
    if z == 0:
        return 1.0, z, walk
    n = len(bits)
    # Practical normal approximation for the NIST cusum statistic.
    p = math.erfc(z / np.sqrt(2 * n))
    return float(p), z, walk


def _pattern_counts_circular(bits, m):
    n = len(bits)
    ext = np.concatenate([bits, bits[:m - 1]])
    counts = np.zeros(2 ** m, dtype=float)
    for i in range(n):
        value = 0
        for b in ext[i:i + m]:
            value = (value << 1) | int(b)
        counts[value] += 1
    return counts


def nist_approximate_entropy(bits, m=2):
    n = len(bits)

    def phi(mm):
        counts = _pattern_counts_circular(bits, mm)
        probs = counts[counts > 0] / n
        return float(np.sum(probs * np.log(probs)))

    ap_en = phi(m) - phi(m + 1)
    chi2 = 2 * n * (np.log(2) - ap_en)
    df = 2 ** (m - 1)
    p = float(stats.chi2.sf(chi2, df))
    return p, float(ap_en), float(chi2), df


def rolling_hurst_rs(series, window=768, step=128):
    """Rolling R/S Hurst procena kroz vreme."""
    x = np.asarray(series, dtype=float)
    centers = []
    hvals = []
    r2vals = []
    for start in range(0, len(x) - window + 1, step):
        seg = x[start:start + window]
        h, _, r2, _, _ = hurst_rs(seg, min_window=8, max_window=max(32, window // 4))
        centers.append(start + window // 2 + 1)
        hvals.append(h)
        r2vals.append(r2)
    return (
        np.asarray(centers, dtype=float),
        np.asarray(hvals, dtype=float),
        np.asarray(r2vals, dtype=float),
    )


def rolling_higuchi_fd(series, window=768, step=128, kmax=32):
    """Rolling Higuchi FD procena kroz vreme."""
    x = np.asarray(series, dtype=float)
    centers = []
    fdvals = []
    r2vals = []
    for start in range(0, len(x) - window + 1, step):
        seg = normalize01(x[start:start + window])
        fd, _, r2, _, _ = higuchi_fd(seg, kmax=min(kmax, window // 4))
        centers.append(start + window // 2 + 1)
        fdvals.append(fd)
        r2vals.append(r2)
    return (
        np.asarray(centers, dtype=float),
        np.asarray(fdvals, dtype=float),
        np.asarray(r2vals, dtype=float),
    )



# ─────────────────────────────────────────────────────────────────────
# KORAK 1: Weierstrass-ova funkcija nad svih 4624 izvucenih kombinacija
#   f(t) = lex-indeks t-tog izvlacenja u prostoru svih 39C7 kombinacija
# ─────────────────────────────────────────────────────────────────────
draws = read_loto_csv(CSV_DRAWS)
N = len(draws)
lex_idx = np.array([lex_rank_1based(c) for c in draws], dtype=np.float64)
t = np.arange(1, N + 1, dtype=np.float64)

print()
print("KarlWeierstrass v2 — KORAK 1: formiranje krive f(t)")
print(f"  CSV:                  {CSV_DRAWS}")
print(f"  Ucitano izvucenja:    {N}")
print(f"  C(39,7):              {TOTAL_COMBOS:,}")
print(f"  Prvih 10 lex-indeksa: {lex_idx[:10].astype(int).tolist()}")
print(f"  min={int(lex_idx.min()):,}  max={int(lex_idx.max()):,}  "
      f"mean={int(lex_idx.mean()):,}  std={int(lex_idx.std()):,}")
print()

# Vizualizacija krive (full-width radi citljivosti)
fig, ax = plt.subplots(figsize=(16, 5))
ax.plot(t, lex_idx, linewidth=0.6, color="steelblue")
ax.set_title(
    f"KORAK 1: Weierstrass-tipa kriva f(t) = lex-indeks  "
    f"({N} izvlacenja, opseg 1..{TOTAL_COMBOS:,})",
    fontsize=12,
)
ax.set_xlabel("t (redni broj izvlacenja)")
ax.set_ylabel("lex-indeks u 39C7")
ax.set_xlim(1, N)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(PNG_PATH, dpi=150, bbox_inches="tight")
plt.show()

elapsed = time.time() - T0

with open(TXT_PATH, "w", encoding="utf-8") as f:
    f.write("KarlWeierstrass algorithm for Loto 7/39 — v2\n")
    f.write("=" * 60 + "\n\n")
    f.write("KORAK 1: Weierstrass-ova funkcija nad svih izvucenih kombinacija\n\n")
    f.write(f"  CSV izvucenih:        {CSV_DRAWS}\n")
    f.write(f"  Prostor svih 39C7:    {CSV_ALL_COMBOS}\n")
    f.write(f"  Ucitano izvucenja:    {N}\n")
    f.write(f"  C(39,7):              {TOTAL_COMBOS:,}\n")
    f.write(f"  PNG:                  {PNG_PATH}\n\n")

    f.write("Definicija krive:\n")
    f.write("  f(t) = lex-indeks izvucene kombinacije u skupu svih 39C7\n")
    f.write("  t    = 1..N (redni broj izvlacenja)\n\n")

    f.write("Prvih 10 tacaka:\n")
    for i in range(min(10, N)):
        f.write(f"  t={i+1:<4} combo={draws[i]}  lex={int(lex_idx[i])}\n")
    f.write("\n")

    f.write("Statistike krive:\n")
    f.write(f"  min  = {int(lex_idx.min()):,}\n")
    f.write(f"  max  = {int(lex_idx.max()):,}\n")
    f.write(f"  mean = {int(lex_idx.mean()):,}\n")
    f.write(f"  std  = {int(lex_idx.std()):,}\n\n")

    f.write(f"Ukupno vreme: {timedelta(seconds=int(elapsed))} ({elapsed:.1f} s)\n")
    f.write("\nKraj KORAKA 1.\n")

print(f"PNG saved → {PNG_PATH}")
print(f"TXT saved → {TXT_PATH}")
print(f"Ukupno vreme: {timedelta(seconds=int(elapsed))} ({elapsed:.1f} s)")
print()
"""
KarlWeierstrass v2 — KORAK 1: formiranje krive f(t)
  CSV:                  /Users/4c/Desktop/GHQ/data/loto7_4624_k43.csv
  Ucitano izvucenja:    4624
  C(39,7):              15,380,937
  Prvih 10 lex-indeksa: [9774196, 3064646, 14623405, 15178170, 5219568, 10167673, 2903730, 9326447, 1398581, 4751586]
  min=6,906  max=15,380,242  mean=7,838,790  std=4,464,599

PNG saved → /1_KarlWeierstrass_v2_1.png
TXT saved → /1_KarlWeierstrass_v2.txt
Ukupno vreme: 0:00:14 (14.8 s)
"""



# ─────────────────────────────────────────────────────────────────────
# KORAK 2a: Aparat — Brownovo kretanje
#   Cilj: izmeriti koliko se f(t) ponasa kao klasicno Brownovo kretanje.
#   - inkrementi dX(t) = f(t+1) - f(t)
#   - stacionarnost: poredjenje prve i druge polovine niza inkremenata
#   - normalnost inkremenata: Shapiro-Wilk + Anderson-Darling
#   - varijansa scaling: Var(X(t+k)-X(t)) ~ k  (Brown: slope = 1 u log-log)
# ─────────────────────────────────────────────────────────────────────
T0_2A = time.time()

incr = np.diff(lex_idx)
half = len(incr) // 2
mean_full, std_full = float(incr.mean()), float(incr.std())
mean_h1, std_h1 = float(incr[:half].mean()), float(incr[:half].std())
mean_h2, std_h2 = float(incr[half:].mean()), float(incr[half:].std())

# Normalnost
rng_2a = np.random.default_rng(42)
n_for_shap = min(5000, len(incr))
sub = rng_2a.choice(incr, n_for_shap, replace=False)
shap_stat, shap_p = stats.shapiro(sub)
ad_res = stats.anderson(incr, dist="norm")
ad_crit_5 = float(ad_res.critical_values[2])  # 5% nivo

# Varijansa scaling Var(X(t+k) - X(t)) vs k
ks_scale = [1, 2, 5, 10, 20, 50, 100, 200, 500]
variances = []
for k in ks_scale:
    if k < N:
        d = lex_idx[k:] - lex_idx[:-k]
        variances.append((k, float(d.var())))
ks_arr = np.array([x[0] for x in variances], dtype=float)
vs_arr = np.array([x[1] for x in variances], dtype=float)
brown_slope, brown_intercept = np.polyfit(np.log(ks_arr), np.log(vs_arr), 1)

print()
print("KORAK 2a: Aparat — Brownovo kretanje")
print(f"  inkrementi:  N={len(incr)}  mean={mean_full:,.1f}  std={std_full:,.1f}")
print(f"  1. polovina: mean={mean_h1:,.1f}  std={std_h1:,.1f}")
print(f"  2. polovina: mean={mean_h2:,.1f}  std={std_h2:,.1f}")
print(f"  Shapiro-Wilk (n={n_for_shap}):  stat={shap_stat:.4f}  p={shap_p:.4f}")
print(f"     ⇒ {'normalno (p>0.05)' if shap_p > 0.05 else 'NIJE normalno (p<=0.05)'}")
print(f"  Anderson-Darling: stat={ad_res.statistic:.3f}  crit(5%)={ad_crit_5:.3f}")
print(f"     ⇒ {'normalno' if ad_res.statistic < ad_crit_5 else 'NIJE normalno'}")
print(f"  Varijansa scaling slope: {brown_slope:.4f}   (Brown: ≈ 1.0)")
print()

# PNG za KORAK 2a: 2x2 panel
fig2a, ax2a = plt.subplots(2, 2, figsize=(14, 9))
fig2a.suptitle(f"KORAK 2a: Brownovo kretanje nad f(t)  (N={N})",
               fontsize=13, fontweight="bold")

ax2a[0, 0].plot(np.arange(1, len(incr) + 1), incr, linewidth=0.4, color="darkorange")
ax2a[0, 0].axhline(0, color="black", linewidth=0.6)
ax2a[0, 0].set_title("Inkrementi dX(t) = f(t+1) - f(t)")
ax2a[0, 0].set_xlabel("t")
ax2a[0, 0].set_ylabel("dX")

ax2a[0, 1].hist(incr, bins=60, color="coral", alpha=0.85, density=True, edgecolor="white")
xs = np.linspace(incr.min(), incr.max(), 400)
ax2a[0, 1].plot(xs, stats.norm.pdf(xs, mean_full, std_full),
                "k--", linewidth=1.4,
                label=f"N(μ={mean_full:,.0f}, σ={std_full:,.0f})")
ax2a[0, 1].set_title(f"Histogram inkremenata  (Shapiro p={shap_p:.3f})")
ax2a[0, 1].legend(fontsize=8)

ax2a[1, 0].loglog(ks_arr, vs_arr, "o-", color="goldenrod")
fit_x = np.array([ks_arr.min(), ks_arr.max()])
ax2a[1, 0].loglog(fit_x, np.exp(brown_intercept) * fit_x ** brown_slope,
                  "k--", label=f"slope={brown_slope:.3f}  (Brown=1)")
ax2a[1, 0].set_title("Varijansa scaling: Var(X(t+k)-X(t)) vs k")
ax2a[1, 0].set_xlabel("k")
ax2a[1, 0].set_ylabel("Var")
ax2a[1, 0].legend(fontsize=8)

stats.probplot(incr, dist="norm", plot=ax2a[1, 1])
ax2a[1, 1].set_title("QQ plot inkremenata vs Normalna")
ax2a[1, 1].get_lines()[0].set_markersize(2.0)

for r in range(2):
    for c in range(2):
        ax2a[r, c].spines["top"].set_visible(False)
        ax2a[r, c].spines["right"].set_visible(False)

fig2a.tight_layout()
fig2a.savefig(PNG_PATH_2A, dpi=150, bbox_inches="tight")
plt.show()

# Append rezultata KORAKA 2a u isti TXT
with open(TXT_PATH, "a", encoding="utf-8") as f:
    f.write("\n")
    f.write("=" * 60 + "\n")
    f.write("KORAK 2a: Aparat — Brownovo kretanje\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"  PNG:                  {PNG_PATH_2A}\n\n")
    f.write("Inkrementi dX(t) = f(t+1) - f(t):\n")
    f.write(f"  N           = {len(incr)}\n")
    f.write(f"  mean        = {mean_full:,.2f}\n")
    f.write(f"  std         = {std_full:,.2f}\n")
    f.write(f"  1. polovina: mean={mean_h1:,.2f}  std={std_h1:,.2f}\n")
    f.write(f"  2. polovina: mean={mean_h2:,.2f}  std={std_h2:,.2f}\n\n")
    f.write("Normalnost inkremenata:\n")
    f.write(f"  Shapiro-Wilk (n={n_for_shap}):  stat={shap_stat:.4f}  p={shap_p:.4f}\n")
    f.write(f"     {'normalno (p>0.05)' if shap_p > 0.05 else 'NIJE normalno (p<=0.05)'}\n")
    f.write(f"  Anderson-Darling:  stat={ad_res.statistic:.3f}  crit(5%)={ad_crit_5:.3f}\n")
    f.write(f"     {'normalno' if ad_res.statistic < ad_crit_5 else 'NIJE normalno'}\n\n")
    f.write("Varijansa scaling Var(X(t+k)-X(t)) vs k  (Brown: linearno, slope ≈ 1):\n")
    f.write(f"  {'k':<8}{'Var':>20}\n")
    for k, v in variances:
        f.write(f"  {k:<8}{v:>20,.2f}\n")
    f.write(f"  log-log slope = {brown_slope:.4f}   (Brown referenca: 1.000)\n\n")

    elapsed_2a = time.time() - T0_2A
    f.write(f"Vreme KORAKA 2a: {timedelta(seconds=int(elapsed_2a))} ({elapsed_2a:.1f} s)\n")
    f.write("\nKraj KORAKA 2a.\n")

print(f"PNG saved → {PNG_PATH_2A}")
print(f"TXT updated → {TXT_PATH}")
print(f"Vreme KORAKA 2a: {timedelta(seconds=int(time.time()-T0_2A))} "
      f"({time.time()-T0_2A:.1f} s)")
print()
"""
KORAK 2a (Aparat: Brownovo kretanje) u 1_KarlWeierstrass_v2.py. Posebna PNG slika, TXT proširen.

KORAK 2a — Brownovo kretanje (inkrementi, stacionarnost po polovinama, Shapiro/Anderson, scaling Var vs k, log-log slope)
           Crta 2x2 panel u 1_KarlWeierstrass_v2_korak2a.png, dopisuje rezultate u isti 1_KarlWeierstrass_v2.txt (append)


KORAK 2a: Aparat — Brownovo kretanje
  inkrementi:  N=4623  mean=-2,003.3  std=6,276,910.0
  1. polovina: mean=-2,315.1  std=6,172,737.3
  2. polovina: mean=-1,691.6  std=6,379,337.8
  Shapiro-Wilk (n=4623):  stat=0.9924  p=0.0000
     ⇒ NIJE normalno (p<=0.05)
  Anderson-Darling: stat=4.702  crit(5%)=0.786
     ⇒ NIJE normalno
  Varijansa scaling slope: 0.0018   (Brown: ≈ 1.0)

PNG saved →   /1_KarlWeierstrass_v2_2a.png
TXT updated → /1_KarlWeierstrass_v2.txt
Vreme KORAKA 2a: 0:00:16 (16.8 s)
"""



# ─────────────────────────────────────────────────────────────────────
# KORAK 2b: Aparat — Hurst eksponent (R/S analiza)
#   Cilj: izmeriti dugorocnu memoriju krive f(t).
#   H ~ 0.5  -> Brown/random referenca
#   H > 0.5  -> perzistentnost
#   H < 0.5  -> anti-perzistentnost
# ─────────────────────────────────────────────────────────────────────
T0_2B = time.time()

hurst_f, hurst_intercept, hurst_r2, hurst_windows, hurst_rs_values = hurst_rs(lex_idx)
hurst_incr, _, hurst_incr_r2, hurst_incr_windows, hurst_incr_rs_values = hurst_rs(incr)

if hurst_f < 0.45:
    hurst_note = "anti-perzistentno (ispod random/Brown reference)"
elif hurst_f > 0.55:
    hurst_note = "perzistentno (iznad random/Brown reference)"
else:
    hurst_note = "blizu random/Brown reference"

print()
print("KORAK 2b: Aparat — Hurst eksponent (R/S analiza)")
print(f"  H(f(t))  = {hurst_f:.4f}   R²={hurst_r2:.4f}   ⇒ {hurst_note}")
print(f"  H(dX(t)) = {hurst_incr:.4f}   R²={hurst_incr_r2:.4f}   (kontrola nad inkrementima)")
print(f"  windows: {hurst_windows.astype(int).tolist()}")
print()

fig2b, ax2b = plt.subplots(1, 2, figsize=(14, 5))
fig2b.suptitle(f"KORAK 2b: Hurst eksponent / R-S analiza  (N={N})",
               fontsize=13, fontweight="bold")

ax2b[0].loglog(hurst_windows, hurst_rs_values, "o-", color="darkslateblue",
               label="R/S f(t)")
fit_y = np.exp(hurst_intercept) * hurst_windows ** hurst_f
ax2b[0].loglog(hurst_windows, fit_y, "k--",
               label=f"H={hurst_f:.3f}, R²={hurst_r2:.3f}")
ax2b[0].set_title("R/S analiza nad f(t)")
ax2b[0].set_xlabel("window")
ax2b[0].set_ylabel("mean R/S")
ax2b[0].legend(fontsize=8)
ax2b[0].grid(True, alpha=0.25, which="both")

ax2b[1].loglog(hurst_incr_windows, hurst_incr_rs_values, "o-", color="seagreen",
               label="R/S dX(t)")
fit_incr = np.exp(np.polyfit(np.log(hurst_incr_windows),
                             np.log(hurst_incr_rs_values), 1)[1]) * (
                                 hurst_incr_windows ** hurst_incr)
ax2b[1].loglog(hurst_incr_windows, fit_incr, "k--",
               label=f"H={hurst_incr:.3f}, R²={hurst_incr_r2:.3f}")
ax2b[1].set_title("Kontrola: R/S analiza nad inkrementima")
ax2b[1].set_xlabel("window")
ax2b[1].set_ylabel("mean R/S")
ax2b[1].legend(fontsize=8)
ax2b[1].grid(True, alpha=0.25, which="both")

for a in ax2b:
    a.spines["top"].set_visible(False)
    a.spines["right"].set_visible(False)

fig2b.tight_layout()
fig2b.savefig(PNG_PATH_2B, dpi=150, bbox_inches="tight")
plt.show()

with open(TXT_PATH, "a", encoding="utf-8") as f:
    f.write("\n")
    f.write("=" * 60 + "\n")
    f.write("KORAK 2b: Aparat — Hurst eksponent (R/S analiza)\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"  PNG:                  {PNG_PATH_2B}\n\n")
    f.write("R/S analiza nad f(t):\n")
    f.write(f"  H(f(t))     = {hurst_f:.4f}\n")
    f.write(f"  R^2         = {hurst_r2:.4f}\n")
    f.write(f"  interpret.  = {hurst_note}\n\n")
    f.write("Kontrola nad inkrementima dX(t):\n")
    f.write(f"  H(dX(t))    = {hurst_incr:.4f}\n")
    f.write(f"  R^2         = {hurst_incr_r2:.4f}\n\n")
    f.write("R/S tacke za f(t):\n")
    f.write(f"  {'window':<10}{'mean R/S':>16}\n")
    for w, rs in zip(hurst_windows.astype(int), hurst_rs_values):
        f.write(f"  {w:<10}{rs:>16,.6f}\n")
    f.write("\n")

    elapsed_2b = time.time() - T0_2B
    f.write(f"Vreme KORAKA 2b: {timedelta(seconds=int(elapsed_2b))} ({elapsed_2b:.1f} s)\n")
    f.write("\nKraj KORAKA 2b.\n")

print(f"PNG saved → {PNG_PATH_2B}")
print(f"TXT updated → {TXT_PATH}")
print(f"Vreme KORAKA 2b: {timedelta(seconds=int(time.time()-T0_2B))} "
      f"({time.time()-T0_2B:.1f} s)")
print()
"""
KORAK 2b (Aparat: Hurst eksponent / R-S analiza) u 1_KarlWeierstrass_v2.py.
Posebna PNG slika, TXT proširen.

KORAK 2b — R/S analiza nad f(t), plus kontrola nad inkrementima dX(t).
           Crta log-log R/S grafikone u 1_KarlWeierstrass_v2_2b.png.


KORAK 2b: Aparat — Hurst eksponent (R/S analiza)
  H(f(t))  = 0.5931   R²=0.9988   ⇒ perzistentno (iznad random/Brown reference)
  H(dX(t)) = 0.0690   R²=0.5986   (kontrola nad inkrementima)
  windows: [8, 12, 18, 27, 40, 59, 86, 125, 182, 264, 383, 556, 807]

PNG saved →   /1_KarlWeierstrass_v2_2b.png
TXT updated → /1_KarlWeierstrass_v2.txt
Vreme KORAKA 2b: 0:00:17 (17.8 s)
"""



# ─────────────────────────────────────────────────────────────────────
# KORAK 2c: Aparat — Fraktalna dimenzija
#   Cilj: izmeriti "hrapavost" krive f(t), bez pretpostavke tangente.
#   Higuchi FD je glavni broj za vremenski niz; Katz FD je kontrolna mera.
# ─────────────────────────────────────────────────────────────────────
T0_2C = time.time()

lex_norm = normalize01(lex_idx)
higuchi_f, higuchi_intercept, higuchi_r2, higuchi_k, higuchi_lk = higuchi_fd(lex_norm)
katz_f = katz_fd(lex_idx)

if higuchi_f < 1.2:
    fd_note = "relativno glatko / slaba fraktalna hrapavost"
elif higuchi_f > 1.7:
    fd_note = "vrlo hrapavo / jako fraktalno-random ponasanje"
else:
    fd_note = "srednja fraktalna hrapavost"

print()
print("KORAK 2c: Aparat — Fraktalna dimenzija")
print(f"  Higuchi FD(f(t)) = {higuchi_f:.4f}   R²={higuchi_r2:.4f}   ⇒ {fd_note}")
print(f"  Katz FD(f(t))    = {katz_f:.4f}   (kontrolna mera)")
print(f"  kmax/ks: {higuchi_k.astype(int).tolist()}")
print()

fig2c, ax2c = plt.subplots(1, 2, figsize=(14, 5))
fig2c.suptitle(f"KORAK 2c: Fraktalna dimenzija nad f(t)  (N={N})",
               fontsize=13, fontweight="bold")

ax2c[0].plot(t, lex_norm, linewidth=0.55, color="steelblue")
ax2c[0].set_title("Normalizovana kriva f(t)")
ax2c[0].set_xlabel("t")
ax2c[0].set_ylabel("normalize01(f(t))")
ax2c[0].grid(True, alpha=0.25)

x_fit = np.log(1.0 / higuchi_k)
y_fit = np.log(higuchi_lk)
ax2c[1].plot(x_fit, y_fit, "o-", color="purple", label="Higuchi L(k)")
ax2c[1].plot(x_fit, higuchi_intercept + higuchi_f * x_fit, "k--",
             label=f"FD={higuchi_f:.3f}, R²={higuchi_r2:.3f}")
ax2c[1].set_title("Higuchi fit: log(L(k)) vs log(1/k)")
ax2c[1].set_xlabel("log(1/k)")
ax2c[1].set_ylabel("log(L(k))")
ax2c[1].legend(fontsize=8)
ax2c[1].grid(True, alpha=0.25)

for a in ax2c:
    a.spines["top"].set_visible(False)
    a.spines["right"].set_visible(False)

fig2c.tight_layout()
fig2c.savefig(PNG_PATH_2C, dpi=150, bbox_inches="tight")
plt.show()

with open(TXT_PATH, "a", encoding="utf-8") as f:
    f.write("\n")
    f.write("=" * 60 + "\n")
    f.write("KORAK 2c: Aparat — Fraktalna dimenzija\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"  PNG:                  {PNG_PATH_2C}\n\n")
    f.write("Fraktalna dimenzija nad f(t):\n")
    f.write(f"  Higuchi FD = {higuchi_f:.4f}\n")
    f.write(f"  R^2        = {higuchi_r2:.4f}\n")
    f.write(f"  Katz FD    = {katz_f:.4f}\n")
    f.write(f"  interpret. = {fd_note}\n\n")
    f.write("Higuchi tacke:\n")
    f.write(f"  {'k':<8}{'L(k)':>18}{'log(1/k)':>18}{'log(L(k))':>18}\n")
    for k, lk in zip(higuchi_k.astype(int), higuchi_lk):
        f.write(f"  {k:<8}{lk:>18,.8f}{np.log(1.0/k):>18,.8f}{np.log(lk):>18,.8f}\n")
    f.write("\n")

    elapsed_2c = time.time() - T0_2C
    f.write(f"Vreme KORAKA 2c: {timedelta(seconds=int(elapsed_2c))} ({elapsed_2c:.1f} s)\n")
    f.write("\nKraj KORAKA 2c.\n")

print(f"PNG saved → {PNG_PATH_2C}")
print(f"TXT updated → {TXT_PATH}")
print(f"Vreme KORAKA 2c: {timedelta(seconds=int(time.time()-T0_2C))} "
      f"({time.time()-T0_2C:.1f} s)")
print()
"""
KORAK 2c (Aparat: Fraktalna dimenzija) u 1_KarlWeierstrass_v2.py.
Posebna PNG slika, TXT proširen.

KORAK 2c — Higuchi FD nad f(t), plus Katz FD kao kontrolna mera.
           Crta normalizovanu krivu i Higuchi log-log fit u 1_KarlWeierstrass_v2_2c.png.


Dodato:
1_KarlWeierstrass_v2_2c.png
Higuchi fractal dimension nad f(t)
Katz FD kao kontrolna mera
TXT append za 2c
komentar-log blok za 2c


KORAK 2c: Aparat — Fraktalna dimenzija
  Higuchi FD(f(t)) = 1.9988   R²=1.0000   ⇒ vrlo hrapavo / jako fraktalno-random ponasanje
  Katz FD(f(t))    = 6.7256   (kontrolna mera)
  kmax/ks: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64]

PNG saved →   /1_KarlWeierstrass_v2_2c.png
TXT updated → /1_KarlWeierstrass_v2.txt
Vreme KORAKA 2c: 0:00:20 (20.6 s)
"""



# ─────────────────────────────────────────────────────────────────────
# KORAK 2a3a: Aparat 2a Brownovo kretanje + Test 3a Hurst eksponent
#   Cilj: proveriti da li Brown-inkrementi imaju Hurst signal iznad
#   shuffled/random Brown reference.
# ─────────────────────────────────────────────────────────────────────
T0_2A3A = time.time()

brown_incr_centered = incr - incr.mean()
brown_path = np.cumsum(brown_incr_centered)

hurst_brown_path, brown_path_intercept, brown_path_r2, brown_path_windows, brown_path_rs = hurst_rs(brown_path)
hurst_brown_incr, brown_incr_intercept, brown_incr_r2, brown_incr_windows, brown_incr_rs = hurst_rs(brown_incr_centered)

rng_2a3a = np.random.default_rng(42)
shuffle_runs = 100
shuffle_h = []
for _ in range(shuffle_runs):
    shuffled = rng_2a3a.permutation(brown_incr_centered)
    shuffled_path = np.cumsum(shuffled)
    h_shuf, _, _, _, _ = hurst_rs(shuffled_path)
    shuffle_h.append(h_shuf)
shuffle_h = np.asarray(shuffle_h, dtype=float)
shuffle_mean = float(shuffle_h.mean())
shuffle_std = float(shuffle_h.std(ddof=1))
shuffle_persist_p = float(np.mean(shuffle_h >= hurst_brown_path))
shuffle_z = (hurst_brown_path - shuffle_mean) / (shuffle_std + 1e-12)

if shuffle_persist_p <= 0.05 and hurst_brown_path > shuffle_mean:
    brown_hurst_note = "Hurst signal iznad shuffled Brown reference"
elif shuffle_persist_p >= 0.95 and hurst_brown_path < shuffle_mean:
    brown_hurst_note = "Hurst signal ispod shuffled Brown reference"
else:
    brown_hurst_note = "nije jak odmak od shuffled Brown reference"

print()
print("KORAK 2a3a: Aparat 2a Brownovo kretanje + Test 3a Hurst eksponent")
print(f"  H(Brown-putanja iz centriranih dX) = {hurst_brown_path:.4f}   R²={brown_path_r2:.4f}")
print(f"  H(centrirani dX)                  = {hurst_brown_incr:.4f}   R²={brown_incr_r2:.4f}")
print(f"  shuffled referenca: mean={shuffle_mean:.4f}  std={shuffle_std:.4f}  z={shuffle_z:.2f}")
print(f"  p(H_shuffled >= H_observed) = {shuffle_persist_p:.4f}   ⇒ {brown_hurst_note}")
print()

fig2a3a, ax2a3a = plt.subplots(1, 3, figsize=(16, 5))
fig2a3a.suptitle("KORAK 2a3a: Brownovo kretanje + Hurst test",
                 fontsize=13, fontweight="bold")

ax2a3a[0].plot(np.arange(1, len(brown_path) + 1), brown_path,
               linewidth=0.55, color="darkorange")
ax2a3a[0].axhline(0, color="black", linewidth=0.6)
ax2a3a[0].set_title("Brown-putanja iz centriranih inkremenata")
ax2a3a[0].set_xlabel("t")
ax2a3a[0].set_ylabel("cumsum(dX - mean(dX))")
ax2a3a[0].grid(True, alpha=0.25)

ax2a3a[1].loglog(brown_path_windows, brown_path_rs, "o-", color="darkslateblue",
                 label="R/S Brown-putanja")
fit_path = np.exp(brown_path_intercept) * brown_path_windows ** hurst_brown_path
ax2a3a[1].loglog(brown_path_windows, fit_path, "k--",
                 label=f"H={hurst_brown_path:.3f}, R²={brown_path_r2:.3f}")
ax2a3a[1].set_title("R/S nad Brown-putanjom")
ax2a3a[1].set_xlabel("window")
ax2a3a[1].set_ylabel("mean R/S")
ax2a3a[1].legend(fontsize=8)
ax2a3a[1].grid(True, alpha=0.25, which="both")

ax2a3a[2].hist(shuffle_h, bins=20, color="lightgray", edgecolor="white")
ax2a3a[2].axvline(hurst_brown_path, color="crimson", linewidth=2,
                  label=f"observed H={hurst_brown_path:.3f}")
ax2a3a[2].axvline(shuffle_mean, color="black", linestyle="--",
                  label=f"shuffle mean={shuffle_mean:.3f}")
ax2a3a[2].set_title("Shuffled Brown H referenca")
ax2a3a[2].set_xlabel("H shuffled path")
ax2a3a[2].set_ylabel("broj")
ax2a3a[2].legend(fontsize=8)

for a in ax2a3a:
    a.spines["top"].set_visible(False)
    a.spines["right"].set_visible(False)

fig2a3a.tight_layout()
fig2a3a.savefig(PNG_PATH_2A3A, dpi=150, bbox_inches="tight")
plt.show()

with open(TXT_PATH, "a", encoding="utf-8") as f:
    f.write("\n")
    f.write("=" * 60 + "\n")
    f.write("KORAK 2a3a: Aparat 2a Brownovo kretanje + Test 3a Hurst eksponent\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"  PNG:                  {PNG_PATH_2A3A}\n\n")
    f.write("Hurst nad Brown aparatom:\n")
    f.write(f"  H(Brown-putanja)      = {hurst_brown_path:.4f}\n")
    f.write(f"  R^2 Brown-putanja     = {brown_path_r2:.4f}\n")
    f.write(f"  H(centrirani dX)      = {hurst_brown_incr:.4f}\n")
    f.write(f"  R^2 centrirani dX     = {brown_incr_r2:.4f}\n\n")
    f.write("Shuffled Brown referenca:\n")
    f.write(f"  runs                  = {shuffle_runs}\n")
    f.write(f"  mean H                = {shuffle_mean:.4f}\n")
    f.write(f"  std H                 = {shuffle_std:.4f}\n")
    f.write(f"  z                     = {shuffle_z:.4f}\n")
    f.write(f"  p(H_shuffled >= H_obs)= {shuffle_persist_p:.4f}\n")
    f.write(f"  interpret.            = {brown_hurst_note}\n\n")
    f.write("R/S tacke za Brown-putanju:\n")
    f.write(f"  {'window':<10}{'mean R/S':>16}\n")
    for w, rs in zip(brown_path_windows.astype(int), brown_path_rs):
        f.write(f"  {w:<10}{rs:>16,.6f}\n")
    f.write("\n")

    elapsed_2a3a = time.time() - T0_2A3A
    f.write(f"Vreme KORAKA 2a3a: {timedelta(seconds=int(elapsed_2a3a))} ({elapsed_2a3a:.1f} s)\n")
    f.write("\nKraj KORAKA 2a3a.\n")

print(f"PNG saved → {PNG_PATH_2A3A}")
print(f"TXT updated → {TXT_PATH}")
print(f"Vreme KORAKA 2a3a: {timedelta(seconds=int(time.time()-T0_2A3A))} "
      f"({time.time()-T0_2A3A:.1f} s)")
print()
"""
KORAK 2a3a (Aparat 2a: Brownovo kretanje + Test 3a: Hurst eksponent)
u 1_KarlWeierstrass_v2.py. Posebna PNG slika, TXT proširen.

KORAK 2a3a — Hurst nad Brown-putanjom iz centriranih inkremenata,
             Hurst nad centriranim inkrementima, plus shuffled Brown referenca.
             Crta rezultat u 1_KarlWeierstrass_v2_2a3a.png.


Dodato:
1_KarlWeierstrass_v2_2a3a.png
Hurst nad Brown-putanjom iz centriranih inkremenata
Hurst nad centriranim inkrementima
shuffled Brown referenca (100 permutacija)
TXT append i komentar-log blok


KORAK 2a3a: Aparat 2a Brownovo kretanje + Test 3a Hurst eksponent
  H(Brown-putanja iz centriranih dX) = 0.6294   R²=0.9950
  H(centrirani dX)                  = 0.0690   R²=0.5986
  shuffled referenca: mean=1.0047  std=0.0098  z=-38.26
  p(H_shuffled >= H_observed) = 1.0000   ⇒ Hurst signal ispod shuffled Brown reference

PNG saved →   /1_KarlWeierstrass_v2_2a3a.png
TXT updated → /1_KarlWeierstrass_v2.txt
Vreme KORAKA 2a3a: 0:00:30 (30.7 s)
"""



# ─────────────────────────────────────────────────────────────────────
# KORAK 2a3b: Aparat 2a Brownovo kretanje + Test 3b Autokorelacija (ACF)
#   Cilj: proveriti da li Brown-inkrementi imaju linearno pamcenje.
#   Za idealni Brown model, inkrementi treba da budu priblizno nekorelisani.
# ─────────────────────────────────────────────────────────────────────
T0_2A3B = time.time()

acf_max_lag = 60
acf_lags = np.arange(0, acf_max_lag + 1)
acf_brown_incr = autocorr_values(brown_incr_centered, acf_max_lag)
acf_brown_path = autocorr_values(brown_path, acf_max_lag)
acf_band = 1.96 / np.sqrt(len(brown_incr_centered))

acf_incr_lags = acf_brown_incr[1:]
max_abs_acf_incr = float(np.max(np.abs(acf_incr_lags)))
sig_lag_count = int(np.sum(np.abs(acf_incr_lags) > acf_band))
top_acf_idx = np.argsort(np.abs(acf_incr_lags))[-10:][::-1] + 1
top_acf_pairs = [(int(lag), float(acf_brown_incr[lag])) for lag in top_acf_idx]

lb_h = 20
lb_q, lb_p = ljung_box_approx(acf_brown_incr, len(brown_incr_centered), lb_h)

rng_2a3b = np.random.default_rng(43)
acf_shuffle_runs = 200
shuffle_max_abs_acf = []
for _ in range(acf_shuffle_runs):
    shuffled = rng_2a3b.permutation(brown_incr_centered)
    shuffled_acf = autocorr_values(shuffled, acf_max_lag)
    shuffle_max_abs_acf.append(float(np.max(np.abs(shuffled_acf[1:]))))
shuffle_max_abs_acf = np.asarray(shuffle_max_abs_acf, dtype=float)
shuffle_acf_mean = float(shuffle_max_abs_acf.mean())
shuffle_acf_std = float(shuffle_max_abs_acf.std(ddof=1))
shuffle_acf_p = float(np.mean(shuffle_max_abs_acf >= max_abs_acf_incr))
shuffle_acf_z = (max_abs_acf_incr - shuffle_acf_mean) / (shuffle_acf_std + 1e-12)

if lb_p <= 0.05 or shuffle_acf_p <= 0.05:
    acf_note = "postoji ACF signal iznad Brown/white-noise reference"
else:
    acf_note = "nema jak ACF signal iznad Brown/white-noise reference"

print()
print("KORAK 2a3b: Aparat 2a Brownovo kretanje + Test 3b Autokorelacija (ACF)")
print(f"  max |ACF(dX)| lag 1..{acf_max_lag}: {max_abs_acf_incr:.4f}")
print(f"  95% band: +/-{acf_band:.4f}   znacajnih lagova: {sig_lag_count}/{acf_max_lag}")
print(f"  Ljung-Box aproks. h={lb_h}: Q={lb_q:.2f}  p={lb_p:.4f}")
print(f"  shuffled max|ACF|: mean={shuffle_acf_mean:.4f} std={shuffle_acf_std:.4f} "
      f"z={shuffle_acf_z:.2f} p={shuffle_acf_p:.4f}")
print(f"  ⇒ {acf_note}")
print()

fig2a3b, ax2a3b = plt.subplots(1, 3, figsize=(16, 5))
fig2a3b.suptitle("KORAK 2a3b: Brownovo kretanje + ACF test",
                 fontsize=13, fontweight="bold")

ax2a3b[0].bar(acf_lags[1:], acf_brown_incr[1:], width=0.8, color="darkorange")
ax2a3b[0].axhline(acf_band, color="crimson", linestyle="--", linewidth=1.2)
ax2a3b[0].axhline(-acf_band, color="crimson", linestyle="--", linewidth=1.2)
ax2a3b[0].axhline(0, color="black", linewidth=0.6)
ax2a3b[0].set_title("ACF centriranih Brown inkremenata")
ax2a3b[0].set_xlabel("lag")
ax2a3b[0].set_ylabel("ACF")

ax2a3b[1].bar(acf_lags[1:], acf_brown_path[1:], width=0.8, color="steelblue")
ax2a3b[1].axhline(acf_band, color="crimson", linestyle="--", linewidth=1.2)
ax2a3b[1].axhline(-acf_band, color="crimson", linestyle="--", linewidth=1.2)
ax2a3b[1].axhline(0, color="black", linewidth=0.6)
ax2a3b[1].set_title("Kontrola: ACF Brown-putanje")
ax2a3b[1].set_xlabel("lag")
ax2a3b[1].set_ylabel("ACF")

ax2a3b[2].hist(shuffle_max_abs_acf, bins=24, color="lightgray", edgecolor="white")
ax2a3b[2].axvline(max_abs_acf_incr, color="crimson", linewidth=2,
                  label=f"observed={max_abs_acf_incr:.3f}")
ax2a3b[2].axvline(shuffle_acf_mean, color="black", linestyle="--",
                  label=f"shuffle mean={shuffle_acf_mean:.3f}")
ax2a3b[2].set_title("Shuffled max |ACF| referenca")
ax2a3b[2].set_xlabel("max |ACF|")
ax2a3b[2].set_ylabel("broj")
ax2a3b[2].legend(fontsize=8)

for a in ax2a3b:
    a.spines["top"].set_visible(False)
    a.spines["right"].set_visible(False)
    a.grid(True, alpha=0.2)

fig2a3b.tight_layout()
fig2a3b.savefig(PNG_PATH_2A3B, dpi=150, bbox_inches="tight")
plt.show()

with open(TXT_PATH, "a", encoding="utf-8") as f:
    f.write("\n")
    f.write("=" * 60 + "\n")
    f.write("KORAK 2a3b: Aparat 2a Brownovo kretanje + Test 3b Autokorelacija (ACF)\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"  PNG:                  {PNG_PATH_2A3B}\n\n")
    f.write("ACF nad centriranim Brown inkrementima:\n")
    f.write(f"  max lag               = {acf_max_lag}\n")
    f.write(f"  95% band              = +/-{acf_band:.6f}\n")
    f.write(f"  max |ACF|             = {max_abs_acf_incr:.6f}\n")
    f.write(f"  znacajnih lagova      = {sig_lag_count}/{acf_max_lag}\n")
    f.write(f"  Ljung-Box h           = {lb_h}\n")
    f.write(f"  Ljung-Box Q           = {lb_q:.6f}\n")
    f.write(f"  Ljung-Box p           = {lb_p:.6f}\n\n")
    f.write("Shuffled max |ACF| referenca:\n")
    f.write(f"  runs                  = {acf_shuffle_runs}\n")
    f.write(f"  mean                  = {shuffle_acf_mean:.6f}\n")
    f.write(f"  std                   = {shuffle_acf_std:.6f}\n")
    f.write(f"  z                     = {shuffle_acf_z:.6f}\n")
    f.write(f"  p(shuffled >= obs)    = {shuffle_acf_p:.6f}\n")
    f.write(f"  interpret.            = {acf_note}\n\n")
    f.write("Top 10 ACF lagova po apsolutnoj vrednosti:\n")
    f.write(f"  {'lag':<8}{'ACF':>16}\n")
    for lag, val in top_acf_pairs:
        f.write(f"  {lag:<8}{val:>16,.8f}\n")
    f.write("\n")

    elapsed_2a3b = time.time() - T0_2A3B
    f.write(f"Vreme KORAKA 2a3b: {timedelta(seconds=int(elapsed_2a3b))} ({elapsed_2a3b:.1f} s)\n")
    f.write("\nKraj KORAKA 2a3b.\n")

print(f"PNG saved → {PNG_PATH_2A3B}")
print(f"TXT updated → {TXT_PATH}")
print(f"Vreme KORAKA 2a3b: {timedelta(seconds=int(time.time()-T0_2A3B))} "
      f"({time.time()-T0_2A3B:.1f} s)")
print()
"""
KORAK 2a3b (Aparat 2a: Brownovo kretanje + Test 3b: Autokorelacija / ACF)
u 1_KarlWeierstrass_v2.py. Posebna PNG slika, TXT proširen.

KORAK 2a3b — ACF nad centriranim Brown inkrementima,
             ACF nad Brown-putanjom kao kontrola,
             Ljung-Box aproksimacija i shuffled max|ACF| referenca.
             Crta rezultat u 1_KarlWeierstrass_v2_2a3b.png.


Aparat 2a Brownovo kretanje + Test 3b Autokorelacija (ACF). Ubacujem novi blok na kraj posle 2a3a, sa PNG 1_KarlWeierstrass_v2_2a3b.png i TXT append,

ACF helper i 2a3b blok. Test je fokusiran na Brown inkremente jer Brown model traži da inkrementi budu što bliže belom šumu; Brown-putanju prikazujem kao kontrolu.

Dodato:
1_KarlWeierstrass_v2_2a3b.png
ACF nad centriranim Brown inkrementima
ACF nad Brown-putanjom kao kontrola
Ljung-Box aproksimacija
shuffled max |ACF| referenca (200 permutacija)
TXT append i komentar-log blok na dnu


KORAK 2a3b: Aparat 2a Brownovo kretanje + Test 3b Autokorelacija (ACF)
  max |ACF(dX)| lag 1..60: 0.5019
  95% band: +/-0.0288   znacajnih lagova: 9/60
  Ljung-Box aproks. h=20: Q=1187.84  p=0.0000
  shuffled max|ACF|: mean=0.0373 std=0.0052 z=88.77 p=0.0000
  ⇒ postoji ACF signal iznad Brown/white-noise reference

PNG saved →   /1_KarlWeierstrass_v2_2a3b.png
TXT updated → /1_KarlWeierstrass_v2.txt
Vreme KORAKA 2a3b: 0:00:15 (15.7 s)
"""



# ─────────────────────────────────────────────────────────────────────
# KORAK 2a3c: Aparat 2a Brownovo kretanje + Test 3c Mutual Information
#   Cilj: proveriti nelinearnu zavisnost Brown-inkremenata po lagovima.
#   ACF hvata linearnu vezu; MI hvata i nelinearnu vezu.
# ─────────────────────────────────────────────────────────────────────
T0_2A3C = time.time()

mi_max_lag = 60
mi_bins = 16
mi_lags = np.arange(1, mi_max_lag + 1)
mi_brown_incr = mutual_information_lags(brown_incr_centered, mi_max_lag, mi_bins)
mi_brown_path = mutual_information_lags(brown_path, mi_max_lag, mi_bins)

max_mi_incr = float(mi_brown_incr.max())
max_mi_lag = int(mi_lags[int(np.argmax(mi_brown_incr))])
top_mi_idx = np.argsort(mi_brown_incr)[-10:][::-1]
top_mi_pairs = [(int(mi_lags[i]), float(mi_brown_incr[i])) for i in top_mi_idx]

rng_2a3c = np.random.default_rng(44)
mi_shuffle_runs = 200
shuffle_max_mi = []
for _ in range(mi_shuffle_runs):
    shuffled = rng_2a3c.permutation(brown_incr_centered)
    shuffled_mi = mutual_information_lags(shuffled, mi_max_lag, mi_bins)
    shuffle_max_mi.append(float(shuffled_mi.max()))
shuffle_max_mi = np.asarray(shuffle_max_mi, dtype=float)
shuffle_mi_mean = float(shuffle_max_mi.mean())
shuffle_mi_std = float(shuffle_max_mi.std(ddof=1))
shuffle_mi_p = float(np.mean(shuffle_max_mi >= max_mi_incr))
shuffle_mi_z = (max_mi_incr - shuffle_mi_mean) / (shuffle_mi_std + 1e-12)

if shuffle_mi_p <= 0.05:
    mi_note = "postoji MI signal iznad shuffled Brown reference"
else:
    mi_note = "nema jak MI signal iznad shuffled Brown reference"

print()
print("KORAK 2a3c: Aparat 2a Brownovo kretanje + Test 3c Mutual Information")
print(f"  max MI(dX) lag 1..{mi_max_lag}: {max_mi_incr:.6f} bits  (lag={max_mi_lag})")
print(f"  shuffled max MI: mean={shuffle_mi_mean:.6f} std={shuffle_mi_std:.6f} "
      f"z={shuffle_mi_z:.2f} p={shuffle_mi_p:.4f}")
print(f"  ⇒ {mi_note}")
print()

fig2a3c, ax2a3c = plt.subplots(1, 3, figsize=(16, 5))
fig2a3c.suptitle("KORAK 2a3c: Brownovo kretanje + Mutual Information test",
                 fontsize=13, fontweight="bold")

ax2a3c[0].plot(mi_lags, mi_brown_incr, "o-", markersize=3, color="darkorange")
ax2a3c[0].set_title("MI centriranih Brown inkremenata")
ax2a3c[0].set_xlabel("lag")
ax2a3c[0].set_ylabel("MI [bits]")
ax2a3c[0].grid(True, alpha=0.25)

ax2a3c[1].plot(mi_lags, mi_brown_path, "o-", markersize=3, color="steelblue")
ax2a3c[1].set_title("Kontrola: MI Brown-putanje")
ax2a3c[1].set_xlabel("lag")
ax2a3c[1].set_ylabel("MI [bits]")
ax2a3c[1].grid(True, alpha=0.25)

ax2a3c[2].hist(shuffle_max_mi, bins=24, color="lightgray", edgecolor="white")
ax2a3c[2].axvline(max_mi_incr, color="crimson", linewidth=2,
                  label=f"observed={max_mi_incr:.4f}")
ax2a3c[2].axvline(shuffle_mi_mean, color="black", linestyle="--",
                  label=f"shuffle mean={shuffle_mi_mean:.4f}")
ax2a3c[2].set_title("Shuffled max MI referenca")
ax2a3c[2].set_xlabel("max MI [bits]")
ax2a3c[2].set_ylabel("broj")
ax2a3c[2].legend(fontsize=8)

for a in ax2a3c:
    a.spines["top"].set_visible(False)
    a.spines["right"].set_visible(False)

fig2a3c.tight_layout()
fig2a3c.savefig(PNG_PATH_2A3C, dpi=150, bbox_inches="tight")
plt.show()

with open(TXT_PATH, "a", encoding="utf-8") as f:
    f.write("\n")
    f.write("=" * 60 + "\n")
    f.write("KORAK 2a3c: Aparat 2a Brownovo kretanje + Test 3c Mutual Information\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"  PNG:                  {PNG_PATH_2A3C}\n\n")
    f.write("Mutual Information nad centriranim Brown inkrementima:\n")
    f.write(f"  max lag               = {mi_max_lag}\n")
    f.write(f"  bins                  = {mi_bins}\n")
    f.write(f"  max MI                = {max_mi_incr:.8f} bits\n")
    f.write(f"  max MI lag            = {max_mi_lag}\n\n")
    f.write("Shuffled max MI referenca:\n")
    f.write(f"  runs                  = {mi_shuffle_runs}\n")
    f.write(f"  mean                  = {shuffle_mi_mean:.8f}\n")
    f.write(f"  std                   = {shuffle_mi_std:.8f}\n")
    f.write(f"  z                     = {shuffle_mi_z:.8f}\n")
    f.write(f"  p(shuffled >= obs)    = {shuffle_mi_p:.8f}\n")
    f.write(f"  interpret.            = {mi_note}\n\n")
    f.write("Top 10 MI lagova:\n")
    f.write(f"  {'lag':<8}{'MI [bits]':>16}\n")
    for lag, val in top_mi_pairs:
        f.write(f"  {lag:<8}{val:>16,.8f}\n")
    f.write("\n")

    elapsed_2a3c = time.time() - T0_2A3C
    f.write(f"Vreme KORAKA 2a3c: {timedelta(seconds=int(elapsed_2a3c))} ({elapsed_2a3c:.1f} s)\n")
    f.write("\nKraj KORAKA 2a3c.\n")

print(f"PNG saved → {PNG_PATH_2A3C}")
print(f"TXT updated → {TXT_PATH}")
print(f"Vreme KORAKA 2a3c: {timedelta(seconds=int(time.time()-T0_2A3C))} "
      f"({time.time()-T0_2A3C:.1f} s)")
print()
"""
KORAK 2a3c (Aparat 2a: Brownovo kretanje + Test 3c: Mutual Information)
u 1_KarlWeierstrass_v2.py. Posebna PNG slika, TXT proširen.

KORAK 2a3c — MI nad centriranim Brown inkrementima,
             MI nad Brown-putanjom kao kontrola,
             shuffled max MI referenca.
             Crta rezultat u 1_KarlWeierstrass_v2_2a3c.png.


Aparat 2a Brownovo kretanje + Test 3c Mutual Information. Fokus je MI nad Brown inkrementima po lagovima, plus Brown-putanja kao kontrola i shuffled referenca; novi PNG _2a3c.

helper za diskretnu Mutual Information (kvantil-binovi, u bitovima) i blok 2a3c. Shuffled referenca služi da odvoji realan nelinearni signal od bias-a zbog binovanja.

Dodato:
1_KarlWeierstrass_v2_2a3c.png
Mutual Information nad centriranim Brown inkrementima
MI nad Brown-putanjom kao kontrola
shuffled max MI referenca (200 permutacija)
TXT append i komentar-log blok


KORAK 2a3c: Aparat 2a Brownovo kretanje + Test 3c Mutual Information
  max MI(dX) lag 1..60: 0.263697 bits  (lag=1)
  shuffled max MI: mean=0.044002 std=0.001862 z=117.97 p=0.0000
  ⇒ postoji MI signal iznad shuffled Brown reference

PNG saved →   /1_KarlWeierstrass_v2_2a3c.png
TXT updated → /1_KarlWeierstrass_v2.txt
Vreme KORAKA 2a3c: 0:00:16 (16.5 s)
"""



# ─────────────────────────────────────────────────────────────────────
# KORAK 2a3d: Aparat 2a Brownovo kretanje + Test 3d Entropy
#   Cilj: izmeriti kompleksnost/nepredvidljivost Brown-inkremenata.
#   Sample entropy radi nad poduzorkom zbog O(n^2); permutation entropy nad celim nizom.
# ─────────────────────────────────────────────────────────────────────
T0_2A3D = time.time()

sampen_m = 2
sampen_r = 0.2
sampen_incr, sampen_a, sampen_b, sampen_n = sample_entropy(
    brown_incr_centered, m=sampen_m, r=sampen_r, max_points=1200
)
sampen_path, _, _, _ = sample_entropy(brown_path, m=sampen_m, r=sampen_r, max_points=1200)

pe_orders = [3, 4, 5]
pe_incr_rows = []
pe_path_rows = []
for order in pe_orders:
    pe, pe_norm, patterns = permutation_entropy(brown_incr_centered, order=order, delay=1)
    pe_incr_rows.append((order, pe, pe_norm, patterns, math.factorial(order)))
    pe_p, pe_p_norm, p_patterns = permutation_entropy(brown_path, order=order, delay=1)
    pe_path_rows.append((order, pe_p, pe_p_norm, p_patterns, math.factorial(order)))

rng_2a3d = np.random.default_rng(45)
entropy_shuffle_runs = 100
shuffle_sampen = []
shuffle_pe4 = []
for _ in range(entropy_shuffle_runs):
    shuffled = rng_2a3d.permutation(brown_incr_centered)
    se, _, _, _ = sample_entropy(shuffled, m=sampen_m, r=sampen_r, max_points=1200)
    _, pe4_norm, _ = permutation_entropy(shuffled, order=4, delay=1)
    if np.isfinite(se):
        shuffle_sampen.append(se)
    shuffle_pe4.append(pe4_norm)
shuffle_sampen = np.asarray(shuffle_sampen, dtype=float)
shuffle_pe4 = np.asarray(shuffle_pe4, dtype=float)

shuffle_sampen_mean = float(shuffle_sampen.mean())
shuffle_sampen_std = float(shuffle_sampen.std(ddof=1))
shuffle_sampen_p_low = float(np.mean(shuffle_sampen <= sampen_incr))
shuffle_sampen_z = (sampen_incr - shuffle_sampen_mean) / (shuffle_sampen_std + 1e-12)

pe4_norm = pe_incr_rows[1][2]
shuffle_pe4_mean = float(shuffle_pe4.mean())
shuffle_pe4_std = float(shuffle_pe4.std(ddof=1))
shuffle_pe4_p_low = float(np.mean(shuffle_pe4 <= pe4_norm))
shuffle_pe4_z = (pe4_norm - shuffle_pe4_mean) / (shuffle_pe4_std + 1e-12)

if shuffle_sampen_p_low <= 0.05 or shuffle_pe4_p_low <= 0.05:
    entropy_note = "entropija je niza od shuffled Brown reference (moguca struktura)"
else:
    entropy_note = "entropija je blizu shuffled Brown reference"

print()
print("KORAK 2a3d: Aparat 2a Brownovo kretanje + Test 3d Sample / Permutation entropy")
print(f"  Sample entropy dX: {sampen_incr:.4f}  (m={sampen_m}, r={sampen_r}, n={sampen_n})")
print(f"  Sample entropy Brown-putanja: {sampen_path:.4f}")
print(f"  Permutation entropy dX order=4: {pe4_norm:.4f} normalizovano")
print(f"  shuffled SampEn: mean={shuffle_sampen_mean:.4f} std={shuffle_sampen_std:.4f} "
      f"z={shuffle_sampen_z:.2f} p_low={shuffle_sampen_p_low:.4f}")
print(f"  shuffled PE4: mean={shuffle_pe4_mean:.4f} std={shuffle_pe4_std:.4f} "
      f"z={shuffle_pe4_z:.2f} p_low={shuffle_pe4_p_low:.4f}")
print(f"  ⇒ {entropy_note}")
print()

fig2a3d, ax2a3d = plt.subplots(1, 3, figsize=(16, 5))
fig2a3d.suptitle("KORAK 2a3d: Brownovo kretanje + Sample / Permutation entropy",
                 fontsize=13, fontweight="bold")

orders = np.array([row[0] for row in pe_incr_rows], dtype=int)
pe_incr_norms = np.array([row[2] for row in pe_incr_rows], dtype=float)
pe_path_norms = np.array([row[2] for row in pe_path_rows], dtype=float)
ax2a3d[0].plot(orders, pe_incr_norms, "o-", color="darkorange", label="dX")
ax2a3d[0].plot(orders, pe_path_norms, "o-", color="steelblue", label="Brown-putanja")
ax2a3d[0].set_ylim(0, 1.05)
ax2a3d[0].set_title("Normalizovana permutation entropy")
ax2a3d[0].set_xlabel("order")
ax2a3d[0].set_ylabel("PE / max PE")
ax2a3d[0].legend(fontsize=8)
ax2a3d[0].grid(True, alpha=0.25)

ax2a3d[1].hist(shuffle_sampen, bins=20, color="lightgray", edgecolor="white")
ax2a3d[1].axvline(sampen_incr, color="crimson", linewidth=2,
                  label=f"observed={sampen_incr:.3f}")
ax2a3d[1].axvline(shuffle_sampen_mean, color="black", linestyle="--",
                  label=f"shuffle mean={shuffle_sampen_mean:.3f}")
ax2a3d[1].set_title("Shuffled Sample Entropy referenca")
ax2a3d[1].set_xlabel("Sample entropy")
ax2a3d[1].set_ylabel("broj")
ax2a3d[1].legend(fontsize=8)

ax2a3d[2].hist(shuffle_pe4, bins=20, color="lightgray", edgecolor="white")
ax2a3d[2].axvline(pe4_norm, color="crimson", linewidth=2,
                  label=f"observed={pe4_norm:.3f}")
ax2a3d[2].axvline(shuffle_pe4_mean, color="black", linestyle="--",
                  label=f"shuffle mean={shuffle_pe4_mean:.3f}")
ax2a3d[2].set_title("Shuffled PE order=4 referenca")
ax2a3d[2].set_xlabel("normalizovana PE")
ax2a3d[2].set_ylabel("broj")
ax2a3d[2].legend(fontsize=8)

for a in ax2a3d:
    a.spines["top"].set_visible(False)
    a.spines["right"].set_visible(False)

fig2a3d.tight_layout()
fig2a3d.savefig(PNG_PATH_2A3D, dpi=150, bbox_inches="tight")
plt.show()

with open(TXT_PATH, "a", encoding="utf-8") as f:
    f.write("\n")
    f.write("=" * 60 + "\n")
    f.write("KORAK 2a3d: Aparat 2a Brownovo kretanje + Test 3d Sample / Permutation entropy\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"  PNG:                  {PNG_PATH_2A3D}\n\n")
    f.write("Sample entropy:\n")
    f.write(f"  m                     = {sampen_m}\n")
    f.write(f"  r                     = {sampen_r}\n")
    f.write(f"  n used                = {sampen_n}\n")
    f.write(f"  SampEn(dX)            = {sampen_incr:.8f}\n")
    f.write(f"  SampEn(Brown-putanja) = {sampen_path:.8f}\n")
    f.write(f"  A count               = {sampen_a}\n")
    f.write(f"  B count               = {sampen_b}\n\n")
    f.write("Permutation entropy:\n")
    f.write(f"  {'order':<8}{'PE bits':>16}{'PE norm':>16}{'patterns':>14}{'max':>10}\n")
    for order, pe, pe_norm, patterns, max_patterns in pe_incr_rows:
        f.write(f"  {order:<8}{pe:>16,.8f}{pe_norm:>16,.8f}{patterns:>14}{max_patterns:>10}\n")
    f.write("\n")
    f.write("Shuffled entropy referenca:\n")
    f.write(f"  runs                  = {entropy_shuffle_runs}\n")
    f.write(f"  SampEn mean           = {shuffle_sampen_mean:.8f}\n")
    f.write(f"  SampEn std            = {shuffle_sampen_std:.8f}\n")
    f.write(f"  SampEn z              = {shuffle_sampen_z:.8f}\n")
    f.write(f"  SampEn p_low          = {shuffle_sampen_p_low:.8f}\n")
    f.write(f"  PE4 mean              = {shuffle_pe4_mean:.8f}\n")
    f.write(f"  PE4 std               = {shuffle_pe4_std:.8f}\n")
    f.write(f"  PE4 z                 = {shuffle_pe4_z:.8f}\n")
    f.write(f"  PE4 p_low             = {shuffle_pe4_p_low:.8f}\n")
    f.write(f"  interpret.            = {entropy_note}\n\n")

    elapsed_2a3d = time.time() - T0_2A3D
    f.write(f"Vreme KORAKA 2a3d: {timedelta(seconds=int(elapsed_2a3d))} ({elapsed_2a3d:.1f} s)\n")
    f.write("\nKraj KORAKA 2a3d.\n")

print(f"PNG saved → {PNG_PATH_2A3D}")
print(f"TXT updated → {TXT_PATH}")
print(f"Vreme KORAKA 2a3d: {timedelta(seconds=int(time.time()-T0_2A3D))} "
      f"({time.time()-T0_2A3D:.1f} s)")
print()
"""
KORAK 2a3d (Aparat 2a: Brownovo kretanje + Test 3d: Sample / Permutation entropy)
u 1_KarlWeierstrass_v2.py. Posebna PNG slika, TXT proširen.

KORAK 2a3d — sample entropy nad centriranim Brown inkrementima,
             permutation entropy nad centriranim Brown inkrementima,
             shuffled entropy referenca.
             Crta rezultat u 1_KarlWeierstrass_v2_2a3d.png.


2a3d: Aparat 2a Brownovo kretanje + Test 3d Sample / Permutation entropy. Dodajem novi blok na kraj, novi PNG _2a3d, TXT append

entropije: sample entropy računam na normalizovanom/smanjenom uzorku radi brzine, permutation entropy na celom nizu; shuffled referenca poredi entropy vrednosti Brown inkremenata sa randomizovanom verzijom.

Dodato:
1_KarlWeierstrass_v2_2a3d.png
Sample entropy nad centriranim Brown inkrementima
Permutation entropy nad centriranim Brown inkrementima
Brown-putanja kao kontrola
shuffled entropy referenca (100 permutacija)
TXT append i komentar-log blok


KORAK 2a3d: Aparat 2a Brownovo kretanje + Test 3d Sample / Permutation entropy
  Sample entropy dX: 2.2117  (m=2, r=0.2, n=1200)
  Sample entropy Brown-putanja: 2.1418
  Permutation entropy dX order=4: 0.9846 normalizovano
  shuffled SampEn: mean=2.2279 std=0.0317 z=-0.51 p_low=0.3100
  shuffled PE4: mean=0.9993 std=0.0002 z=-65.27 p_low=0.0000
  ⇒ entropija je niza od shuffled Brown reference (moguca struktura)

PNG saved →   /1_KarlWeierstrass_v2_2a3d.png
TXT updated → /1_KarlWeierstrass_v2.txt
Vreme KORAKA 2a3d: 0:00:18 (18.6 s)
"""



# ─────────────────────────────────────────────────────────────────────
# KORAK 2a3e: Aparat 2a Brownovo kretanje + Test 3e NIST baterija
#   Cilj: proveriti binarnu nasumicnost Brown-inkremenata.
#   Ovo je prakticna NIST-style baterija, ne kompletan zvanicni NIST paket.
# ─────────────────────────────────────────────────────────────────────
T0_2A3E = time.time()

nist_bits = nist_bits_from_series(brown_incr_centered)
nist_n = len(nist_bits)
nist_ones = int(nist_bits.sum())
nist_zeros = int(nist_n - nist_ones)

monobit_p, monobit_s = nist_monobit(nist_bits)
runs_p, runs_count, runs_pi = nist_runs(nist_bits)
block_p, block_chi2, block_count = nist_block_frequency(nist_bits, block_size=128)
cusum_p, cusum_z, cusum_walk = nist_cumulative_sums(nist_bits)
apen_p, apen_value, apen_chi2, apen_df = nist_approximate_entropy(nist_bits, m=2)

nist_rows = [
    ("Monobit frequency", monobit_p, monobit_s),
    ("Runs", runs_p, runs_count),
    ("Block frequency", block_p, block_chi2),
    ("Cumulative sums", cusum_p, cusum_z),
    ("Approx entropy", apen_p, apen_value),
]
nist_pass_count = int(sum(p > 0.05 for _, p, _ in nist_rows if np.isfinite(p)))
nist_total = int(sum(np.isfinite(p) for _, p, _ in nist_rows))

if nist_pass_count == nist_total:
    nist_note = "svi NIST-style testovi prolaze na 0.05"
elif nist_pass_count >= max(1, nist_total - 1):
    nist_note = "uglavnom prolazi, postoji slab signal za proveru"
else:
    nist_note = "vise NIST-style testova pada (moguća struktura / odstupanje)"

print()
print("KORAK 2a3e: Aparat 2a Brownovo kretanje + Test 3e NIST baterija")
print(f"  bits: n={nist_n}  zeros={nist_zeros}  ones={nist_ones}")
for name, p, stat_val in nist_rows:
    print(f"  {name:<20} p={p:.6f}  stat={stat_val}")
print(f"  prolaz: {nist_pass_count}/{nist_total}  ⇒ {nist_note}")
print()

fig2a3e, ax2a3e = plt.subplots(1, 3, figsize=(16, 5))
fig2a3e.suptitle("KORAK 2a3e: Brownovo kretanje + NIST-style testovi",
                 fontsize=13, fontweight="bold")

ax2a3e[0].bar(["0", "1"], [nist_zeros, nist_ones], color=["steelblue", "darkorange"])
ax2a3e[0].set_title("Binarizovani Brown inkrementi")
ax2a3e[0].set_xlabel("bit")
ax2a3e[0].set_ylabel("broj")
ax2a3e[0].grid(True, alpha=0.2, axis="y")

names = [row[0] for row in nist_rows]
pvals = np.array([row[1] for row in nist_rows], dtype=float)
colors = ["seagreen" if p > 0.05 else "crimson" for p in pvals]
ax2a3e[1].barh(names, pvals, color=colors)
ax2a3e[1].axvline(0.05, color="black", linestyle="--", linewidth=1.2)
ax2a3e[1].set_xlim(0, 1)
ax2a3e[1].set_title("NIST-style p-vrednosti")
ax2a3e[1].set_xlabel("p-value")
ax2a3e[1].grid(True, alpha=0.2, axis="x")

ax2a3e[2].plot(np.arange(1, nist_n + 1), cusum_walk, linewidth=0.7, color="purple")
ax2a3e[2].axhline(0, color="black", linewidth=0.6)
ax2a3e[2].set_title(f"Cumulative sums walk (z={cusum_z})")
ax2a3e[2].set_xlabel("t")
ax2a3e[2].set_ylabel("cum sum")
ax2a3e[2].grid(True, alpha=0.25)

for a in ax2a3e:
    a.spines["top"].set_visible(False)
    a.spines["right"].set_visible(False)

fig2a3e.tight_layout()
fig2a3e.savefig(PNG_PATH_2A3E, dpi=150, bbox_inches="tight")
plt.show()

with open(TXT_PATH, "a", encoding="utf-8") as f:
    f.write("\n")
    f.write("=" * 60 + "\n")
    f.write("KORAK 2a3e: Aparat 2a Brownovo kretanje + Test 3e NIST baterija\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"  PNG:                  {PNG_PATH_2A3E}\n\n")
    f.write("Binarizacija Brown inkremenata:\n")
    f.write("  bit = 1 ako je centrirani dX iznad medijane, inace 0\n")
    f.write(f"  n bits                = {nist_n}\n")
    f.write(f"  zeros                 = {nist_zeros}\n")
    f.write(f"  ones                  = {nist_ones}\n\n")
    f.write("NIST-style testovi (prolaz ako p > 0.05):\n")
    f.write(f"  {'test':<22}{'p-value':>14}{'stat':>18}{'pass':>10}\n")
    for name, p, stat_val in nist_rows:
        f.write(f"  {name:<22}{p:>14,.8f}{float(stat_val):>18,.8f}{str(p > 0.05):>10}\n")
    f.write("\n")
    f.write("Detalji:\n")
    f.write(f"  Runs pi               = {runs_pi:.8f}\n")
    f.write(f"  Block count           = {block_count}\n")
    f.write(f"  Approx entropy m      = 2\n")
    f.write(f"  Approx entropy chi2   = {apen_chi2:.8f}\n")
    f.write(f"  Approx entropy df     = {apen_df}\n")
    f.write(f"  pass count            = {nist_pass_count}/{nist_total}\n")
    f.write(f"  interpret.            = {nist_note}\n\n")

    elapsed_2a3e = time.time() - T0_2A3E
    f.write(f"Vreme KORAKA 2a3e: {timedelta(seconds=int(elapsed_2a3e))} ({elapsed_2a3e:.1f} s)\n")
    f.write("\nKraj KORAKA 2a3e.\n")

print(f"PNG saved → {PNG_PATH_2A3E}")
print(f"TXT updated → {TXT_PATH}")
print(f"Vreme KORAKA 2a3e: {timedelta(seconds=int(time.time()-T0_2A3E))} "
      f"({time.time()-T0_2A3E:.1f} s)")
print()
"""
KORAK 2a3e (Aparat 2a: Brownovo kretanje + Test 3e: NIST baterija)
u 1_KarlWeierstrass_v2.py. Posebna PNG slika, TXT proširen.

KORAK 2a3e — NIST-style testovi nad binarizovanim Brown inkrementima:
             monobit, runs, block frequency, cumulative sums,
             approximate entropy.
             Crta rezultat u 1_KarlWeierstrass_v2_2a3e.png.


2a3e: Aparat 2a Brownovo kretanje + Test 3e NIST baterija nasumičnosti. Dodaću osnovne NIST-style testove nad binarizovanim Brown inkrementima, novi PNG _2a3e, TXT append.

NIST-style helper funkcije pre KORAKA 1, pa zatim blok 2a3e na kraj. Testovi: monobit, runs, block frequency, cumulative sums i approximate entropy.

Dodato:
1_KarlWeierstrass_v2_2a3e.png
NIST-style testovi nad binarizovanim Brown inkrementima:
Monobit frequency
Runs
Block frequency
Cumulative sums
Approximate entropy
TXT append i komentar-log blok


KORAK 2a3e: Aparat 2a Brownovo kretanje + Test 3e NIST baterija
  bits: n=4623  zeros=2312  ones=2311
  Monobit frequency    p=0.988266  stat=0.014707472779848216
  Runs                 p=0.000000  stat=3087
  Block frequency      p=0.999997  stat=9.5
  Cumulative sums      p=0.680479  stat=28
  Approx entropy       p=0.000000  stat=0.6297857663526234
  prolaz: 3/5  ⇒ vise NIST-style testova pada (moguća struktura / odstupanje)

PNG saved →   /1_KarlWeierstrass_v2_2a3e.png
TXT updated → /1_KarlWeierstrass_v2.txt
Vreme KORAKA 2a3e: 0:00:18 (18.6 s)
"""



# ─────────────────────────────────────────────────────────────────────
# KORAK 2b3a: Aparat 2b Hurst/R-S + Test 3a Hurst eksponent
#   Cilj: proveriti da li globalni H(f(t)) iz 2b ostaje stabilan kroz vreme
#   i da li odstupa od shuffled Hurst reference.
# ─────────────────────────────────────────────────────────────────────
T0_2B3A = time.time()

rolling_window = 768
rolling_step = 128
roll_centers, roll_h, roll_r2 = rolling_hurst_rs(
    lex_idx, window=rolling_window, step=rolling_step
)

rng_2b3a = np.random.default_rng(46)
hurst_shuffle_runs = 100
shuffle_h_f = []
for _ in range(hurst_shuffle_runs):
    shuffled_f = rng_2b3a.permutation(lex_idx)
    h_shuf, _, _, _, _ = hurst_rs(shuffled_f)
    shuffle_h_f.append(h_shuf)
shuffle_h_f = np.asarray(shuffle_h_f, dtype=float)

shuffle_h_mean = float(shuffle_h_f.mean())
shuffle_h_std = float(shuffle_h_f.std(ddof=1))
shuffle_h_p_high = float(np.mean(shuffle_h_f >= hurst_f))
shuffle_h_p_low = float(np.mean(shuffle_h_f <= hurst_f))
shuffle_h_z = (hurst_f - shuffle_h_mean) / (shuffle_h_std + 1e-12)

roll_h_mean = float(roll_h.mean())
roll_h_std = float(roll_h.std(ddof=1))
roll_h_min = float(roll_h.min())
roll_h_max = float(roll_h.max())
roll_persistent_count = int(np.sum(roll_h > 0.55))
roll_antipersistent_count = int(np.sum(roll_h < 0.45))

if shuffle_h_p_high <= 0.05 and hurst_f > shuffle_h_mean:
    hurst_2b3a_note = "globalni H je iznad shuffled reference"
elif shuffle_h_p_low <= 0.05 and hurst_f < shuffle_h_mean:
    hurst_2b3a_note = "globalni H je ispod shuffled reference"
else:
    hurst_2b3a_note = "globalni H nije jak odmak od shuffled reference"

print()
print("KORAK 2b3a: Aparat 2b Hurst/R-S + Test 3a Hurst eksponent")
print(f"  global H(f(t)) = {hurst_f:.4f}   R²={hurst_r2:.4f}")
print(f"  rolling H: mean={roll_h_mean:.4f} std={roll_h_std:.4f} "
      f"min={roll_h_min:.4f} max={roll_h_max:.4f}")
print(f"  rolling prozori: persistent={roll_persistent_count}/{len(roll_h)}  "
      f"anti={roll_antipersistent_count}/{len(roll_h)}")
print(f"  shuffled H: mean={shuffle_h_mean:.4f} std={shuffle_h_std:.4f} "
      f"z={shuffle_h_z:.2f} p_high={shuffle_h_p_high:.4f}")
print(f"  ⇒ {hurst_2b3a_note}")
print()

fig2b3a, ax2b3a = plt.subplots(1, 3, figsize=(16, 5))
fig2b3a.suptitle("KORAK 2b3a: Hurst/R-S aparat + Hurst test",
                 fontsize=13, fontweight="bold")

ax2b3a[0].plot(roll_centers, roll_h, "o-", markersize=3, color="darkslateblue")
ax2b3a[0].axhline(0.5, color="black", linestyle="--", linewidth=1.1, label="H=0.5")
ax2b3a[0].axhline(0.55, color="seagreen", linestyle=":", linewidth=1.0)
ax2b3a[0].axhline(0.45, color="crimson", linestyle=":", linewidth=1.0)
ax2b3a[0].set_title(f"Rolling Hurst (window={rolling_window}, step={rolling_step})")
ax2b3a[0].set_xlabel("t centar prozora")
ax2b3a[0].set_ylabel("H")
ax2b3a[0].legend(fontsize=8)
ax2b3a[0].grid(True, alpha=0.25)

ax2b3a[1].hist(shuffle_h_f, bins=22, color="lightgray", edgecolor="white")
ax2b3a[1].axvline(hurst_f, color="crimson", linewidth=2,
                  label=f"observed H={hurst_f:.3f}")
ax2b3a[1].axvline(shuffle_h_mean, color="black", linestyle="--",
                  label=f"shuffle mean={shuffle_h_mean:.3f}")
ax2b3a[1].set_title("Shuffled Hurst referenca")
ax2b3a[1].set_xlabel("H shuffled f(t)")
ax2b3a[1].set_ylabel("broj")
ax2b3a[1].legend(fontsize=8)

ax2b3a[2].plot(roll_centers, roll_r2, "o-", markersize=3, color="steelblue")
ax2b3a[2].set_ylim(0, 1.05)
ax2b3a[2].set_title("Kvalitet rolling R/S fit-a")
ax2b3a[2].set_xlabel("t centar prozora")
ax2b3a[2].set_ylabel("R²")
ax2b3a[2].grid(True, alpha=0.25)

for a in ax2b3a:
    a.spines["top"].set_visible(False)
    a.spines["right"].set_visible(False)

fig2b3a.tight_layout()
fig2b3a.savefig(PNG_PATH_2B3A, dpi=150, bbox_inches="tight")
plt.show()

with open(TXT_PATH, "a", encoding="utf-8") as f:
    f.write("\n")
    f.write("=" * 60 + "\n")
    f.write("KORAK 2b3a: Aparat 2b Hurst/R-S + Test 3a Hurst eksponent\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"  PNG:                  {PNG_PATH_2B3A}\n\n")
    f.write("Globalni Hurst nad f(t):\n")
    f.write(f"  H(f(t))               = {hurst_f:.8f}\n")
    f.write(f"  R^2                   = {hurst_r2:.8f}\n\n")
    f.write("Rolling/local Hurst:\n")
    f.write(f"  window                = {rolling_window}\n")
    f.write(f"  step                  = {rolling_step}\n")
    f.write(f"  broj prozora          = {len(roll_h)}\n")
    f.write(f"  mean H                = {roll_h_mean:.8f}\n")
    f.write(f"  std H                 = {roll_h_std:.8f}\n")
    f.write(f"  min H                 = {roll_h_min:.8f}\n")
    f.write(f"  max H                 = {roll_h_max:.8f}\n")
    f.write(f"  persistent H>0.55     = {roll_persistent_count}/{len(roll_h)}\n")
    f.write(f"  anti H<0.45           = {roll_antipersistent_count}/{len(roll_h)}\n\n")
    f.write("Shuffled Hurst referenca:\n")
    f.write(f"  runs                  = {hurst_shuffle_runs}\n")
    f.write(f"  mean                  = {shuffle_h_mean:.8f}\n")
    f.write(f"  std                   = {shuffle_h_std:.8f}\n")
    f.write(f"  z                     = {shuffle_h_z:.8f}\n")
    f.write(f"  p_high                = {shuffle_h_p_high:.8f}\n")
    f.write(f"  p_low                 = {shuffle_h_p_low:.8f}\n")
    f.write(f"  interpret.            = {hurst_2b3a_note}\n\n")
    f.write("Rolling H tacke:\n")
    f.write(f"  {'center':<10}{'H':>16}{'R^2':>16}\n")
    for center, hval, r2val in zip(roll_centers.astype(int), roll_h, roll_r2):
        f.write(f"  {center:<10}{hval:>16,.8f}{r2val:>16,.8f}\n")
    f.write("\n")

    elapsed_2b3a = time.time() - T0_2B3A
    f.write(f"Vreme KORAKA 2b3a: {timedelta(seconds=int(elapsed_2b3a))} ({elapsed_2b3a:.1f} s)\n")
    f.write("\nKraj KORAKA 2b3a.\n")

print(f"PNG saved → {PNG_PATH_2B3A}")
print(f"TXT updated → {TXT_PATH}")
print(f"Vreme KORAKA 2b3a: {timedelta(seconds=int(time.time()-T0_2B3A))} "
      f"({time.time()-T0_2B3A:.1f} s)")
print()
"""
KORAK 2b3a (Aparat 2b: Hurst/R-S + Test 3a: Hurst eksponent)
u 1_KarlWeierstrass_v2.py. Posebna PNG slika, TXT proširen.

KORAK 2b3a — globalni Hurst nad f(t),
             rolling/local Hurst kroz vreme,
             shuffled Hurst referenca.
             Crta rezultat u 1_KarlWeierstrass_v2_2b3a.png.


2b3a: Aparat 2b Hurst/R-S + Test 3a Hurst eksponent. Ovo radim kao detaljniju Hurst proveru nad f(t): globalni H, rolling/local H po prozorima i shuffled referenca, sa novim PNG _2b3a i TXT append.

2b3a kao Hurst validaciju samog Hurst aparata: globalni H iz 2b, rolling H kroz vreme i shuffled H referenca.

Dodato:
1_KarlWeierstrass_v2_2b3a.png
globalni Hurst nad f(t)
rolling/local Hurst kroz vreme
shuffled Hurst referenca (100 permutacija)
TXT append i komentar-log blok


KORAK 2b3a: Aparat 2b Hurst/R-S + Test 3a Hurst eksponent
  global H(f(t)) = 0.5931   R²=0.9988
  rolling H: mean=0.5974 std=0.0252 min=0.5494 max=0.6466
  rolling prozori: persistent=30/31  anti=0/31
  shuffled H: mean=0.5642 std=0.0173 z=1.67 p_high=0.0300
  ⇒ globalni H je iznad shuffled reference

PNG saved →   /1_KarlWeierstrass_v2_2b3a.png
TXT updated → /1_KarlWeierstrass_v2.txt
Vreme KORAKA 2b3a: 0:00:35 (35.0 s)
"""



# ─────────────────────────────────────────────────────────────────────
# KORAK 2b3b: Aparat 2b Hurst/R-S + Test 3b Autokorelacija (ACF)
#   Cilj: proveriti da li se lokalni Hurst rezimi grupisu kroz vreme.
#   ACF rolling H niza meri memoriju samog Hurst aparata.
# ─────────────────────────────────────────────────────────────────────
T0_2B3B = time.time()

roll_acf_max_lag = min(20, max(1, len(roll_h) - 2))
roll_acf_lags = np.arange(0, roll_acf_max_lag + 1)
acf_roll_h = autocorr_values(roll_h, roll_acf_max_lag)
acf_f_control = autocorr_values(lex_idx, acf_max_lag)

roll_acf_band = 1.96 / np.sqrt(len(roll_h))
roll_acf_body = acf_roll_h[1:]
roll_max_abs_acf = float(np.max(np.abs(roll_acf_body)))
roll_sig_lag_count = int(np.sum(np.abs(roll_acf_body) > roll_acf_band))
roll_top_idx = np.argsort(np.abs(roll_acf_body))[-min(10, len(roll_acf_body)):][::-1] + 1
roll_top_acf_pairs = [(int(lag), float(acf_roll_h[lag])) for lag in roll_top_idx]

roll_lb_h = min(10, roll_acf_max_lag)
roll_lb_q, roll_lb_p = ljung_box_approx(acf_roll_h, len(roll_h), roll_lb_h)

rng_2b3b = np.random.default_rng(47)
roll_acf_shuffle_runs = 500
shuffle_roll_max_abs_acf = []
for _ in range(roll_acf_shuffle_runs):
    shuffled_h = rng_2b3b.permutation(roll_h)
    shuffled_acf = autocorr_values(shuffled_h, roll_acf_max_lag)
    shuffle_roll_max_abs_acf.append(float(np.max(np.abs(shuffled_acf[1:]))))
shuffle_roll_max_abs_acf = np.asarray(shuffle_roll_max_abs_acf, dtype=float)
shuffle_roll_acf_mean = float(shuffle_roll_max_abs_acf.mean())
shuffle_roll_acf_std = float(shuffle_roll_max_abs_acf.std(ddof=1))
shuffle_roll_acf_p = float(np.mean(shuffle_roll_max_abs_acf >= roll_max_abs_acf))
shuffle_roll_acf_z = (
    (roll_max_abs_acf - shuffle_roll_acf_mean) / (shuffle_roll_acf_std + 1e-12)
)

if roll_lb_p <= 0.05 or shuffle_roll_acf_p <= 0.05:
    roll_acf_note = "rolling H ima ACF signal iznad shuffled reference"
else:
    roll_acf_note = "rolling H nema jak ACF signal iznad shuffled reference"

print()
print("KORAK 2b3b: Aparat 2b Hurst/R-S + Test 3b Autokorelacija (ACF)")
print(f"  rolling H max |ACF| lag 1..{roll_acf_max_lag}: {roll_max_abs_acf:.4f}")
print(f"  95% band: +/-{roll_acf_band:.4f}   znacajnih lagova: "
      f"{roll_sig_lag_count}/{roll_acf_max_lag}")
print(f"  Ljung-Box aproks. h={roll_lb_h}: Q={roll_lb_q:.2f}  p={roll_lb_p:.4f}")
print(f"  shuffled max|ACF|: mean={shuffle_roll_acf_mean:.4f} std={shuffle_roll_acf_std:.4f} "
      f"z={shuffle_roll_acf_z:.2f} p={shuffle_roll_acf_p:.4f}")
print(f"  ⇒ {roll_acf_note}")
print()

fig2b3b, ax2b3b = plt.subplots(1, 3, figsize=(16, 5))
fig2b3b.suptitle("KORAK 2b3b: Hurst/R-S aparat + ACF test",
                 fontsize=13, fontweight="bold")

ax2b3b[0].bar(roll_acf_lags[1:], acf_roll_h[1:], width=0.8, color="darkslateblue")
ax2b3b[0].axhline(roll_acf_band, color="crimson", linestyle="--", linewidth=1.2)
ax2b3b[0].axhline(-roll_acf_band, color="crimson", linestyle="--", linewidth=1.2)
ax2b3b[0].axhline(0, color="black", linewidth=0.6)
ax2b3b[0].set_title("ACF rolling/local Hurst niza")
ax2b3b[0].set_xlabel("lag")
ax2b3b[0].set_ylabel("ACF")

ax2b3b[1].bar(np.arange(1, acf_max_lag + 1), acf_f_control[1:], width=0.8,
              color="steelblue")
ax2b3b[1].axhline(1.96 / np.sqrt(N), color="crimson", linestyle="--", linewidth=1.2)
ax2b3b[1].axhline(-1.96 / np.sqrt(N), color="crimson", linestyle="--", linewidth=1.2)
ax2b3b[1].axhline(0, color="black", linewidth=0.6)
ax2b3b[1].set_title("Kontrola: ACF f(t)")
ax2b3b[1].set_xlabel("lag")
ax2b3b[1].set_ylabel("ACF")

ax2b3b[2].hist(shuffle_roll_max_abs_acf, bins=24, color="lightgray", edgecolor="white")
ax2b3b[2].axvline(roll_max_abs_acf, color="crimson", linewidth=2,
                  label=f"observed={roll_max_abs_acf:.3f}")
ax2b3b[2].axvline(shuffle_roll_acf_mean, color="black", linestyle="--",
                  label=f"shuffle mean={shuffle_roll_acf_mean:.3f}")
ax2b3b[2].set_title("Shuffled rolling H max |ACF|")
ax2b3b[2].set_xlabel("max |ACF|")
ax2b3b[2].set_ylabel("broj")
ax2b3b[2].legend(fontsize=8)

for a in ax2b3b:
    a.spines["top"].set_visible(False)
    a.spines["right"].set_visible(False)
    a.grid(True, alpha=0.2)

fig2b3b.tight_layout()
fig2b3b.savefig(PNG_PATH_2B3B, dpi=150, bbox_inches="tight")
plt.show()

with open(TXT_PATH, "a", encoding="utf-8") as f:
    f.write("\n")
    f.write("=" * 60 + "\n")
    f.write("KORAK 2b3b: Aparat 2b Hurst/R-S + Test 3b Autokorelacija (ACF)\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"  PNG:                  {PNG_PATH_2B3B}\n\n")
    f.write("ACF nad rolling/local Hurst nizom:\n")
    f.write(f"  broj rolling H tacaka = {len(roll_h)}\n")
    f.write(f"  max lag               = {roll_acf_max_lag}\n")
    f.write(f"  95% band              = +/-{roll_acf_band:.6f}\n")
    f.write(f"  max |ACF|             = {roll_max_abs_acf:.6f}\n")
    f.write(f"  znacajnih lagova      = {roll_sig_lag_count}/{roll_acf_max_lag}\n")
    f.write(f"  Ljung-Box h           = {roll_lb_h}\n")
    f.write(f"  Ljung-Box Q           = {roll_lb_q:.6f}\n")
    f.write(f"  Ljung-Box p           = {roll_lb_p:.6f}\n\n")
    f.write("Shuffled rolling H max |ACF| referenca:\n")
    f.write(f"  runs                  = {roll_acf_shuffle_runs}\n")
    f.write(f"  mean                  = {shuffle_roll_acf_mean:.6f}\n")
    f.write(f"  std                   = {shuffle_roll_acf_std:.6f}\n")
    f.write(f"  z                     = {shuffle_roll_acf_z:.6f}\n")
    f.write(f"  p(shuffled >= obs)    = {shuffle_roll_acf_p:.6f}\n")
    f.write(f"  interpret.            = {roll_acf_note}\n\n")
    f.write("Top ACF lagovi rolling H po apsolutnoj vrednosti:\n")
    f.write(f"  {'lag':<8}{'ACF':>16}\n")
    for lag, val in roll_top_acf_pairs:
        f.write(f"  {lag:<8}{val:>16,.8f}\n")
    f.write("\n")

    elapsed_2b3b = time.time() - T0_2B3B
    f.write(f"Vreme KORAKA 2b3b: {timedelta(seconds=int(elapsed_2b3b))} ({elapsed_2b3b:.1f} s)\n")
    f.write("\nKraj KORAKA 2b3b.\n")

print(f"PNG saved → {PNG_PATH_2B3B}")
print(f"TXT updated → {TXT_PATH}")
print(f"Vreme KORAKA 2b3b: {timedelta(seconds=int(time.time()-T0_2B3B))} "
      f"({time.time()-T0_2B3B:.1f} s)")
print()
"""
KORAK 2b3b (Aparat 2b: Hurst/R-S + Test 3b: Autokorelacija / ACF)
u 1_KarlWeierstrass_v2.py. Posebna PNG slika, TXT proširen.

KORAK 2b3b — ACF nad rolling/local Hurst nizom,
             ACF nad f(t) kao kontrola,
             Ljung-Box i shuffled max|ACF| referenca.
             Crta rezultat u 1_KarlWeierstrass_v2_2b3b.png.


2b3b: Aparat 2b Hurst/R-S + Test 3b Autokorelacija (ACF). Radiću ACF nad rolling/local Hurst nizom, plus Ljung-Box i shuffled referencu; novi PNG _2b3b, TXT append.

2b3b: ACF rolling H režima, ACF f(t) kao kontrolu, Ljung-Box i shuffled max|ACF| referencu za rolling H.

Dodato:
1_KarlWeierstrass_v2_2b3b.png
ACF nad rolling/local Hurst nizom
ACF nad f(t) kao kontrola
Ljung-Box aproksimacija
shuffled max |ACF| referenca (500 permutacija)
TXT append i komentar-log blok


KORAK 2b3b: Aparat 2b Hurst/R-S + Test 3b Autokorelacija (ACF)
  rolling H max |ACF| lag 1..20: 0.4771
  95% band: +/-0.3520   znacajnih lagova: 1/20
  Ljung-Box aproks. h=10: Q=21.93  p=0.0154
  shuffled max|ACF|: mean=0.3101 std=0.0738 z=2.26 p=0.0280
  ⇒ rolling H ima ACF signal iznad shuffled reference

PNG saved →   /1_KarlWeierstrass_v2_2b3b.png
TXT updated → /1_KarlWeierstrass_v2.txt
Vreme KORAKA 2b3b: 0:00:13 (13.6 s)
"""



# ─────────────────────────────────────────────────────────────────────
# KORAK 2b3c: Aparat 2b Hurst/R-S + Test 3c Mutual Information
#   Cilj: proveriti nelinearnu zavisnost lokalnih Hurst rezima kroz vreme.
#   MI nad rolling H hvata veze koje ACF moze da promasi.
# ─────────────────────────────────────────────────────────────────────
T0_2B3C = time.time()

roll_mi_max_lag = min(12, max(1, len(roll_h) // 2))
roll_mi_bins = min(8, max(3, len(roll_h) // 4))
roll_mi_lags = np.arange(1, roll_mi_max_lag + 1)
mi_roll_h = mutual_information_lags(roll_h, roll_mi_max_lag, roll_mi_bins)

mi_f_control_max_lag = 60
mi_f_control_bins = 16
mi_f_control = mutual_information_lags(lex_idx, mi_f_control_max_lag, mi_f_control_bins)
mi_f_control_lags = np.arange(1, mi_f_control_max_lag + 1)

roll_max_mi = float(mi_roll_h.max())
roll_max_mi_lag = int(roll_mi_lags[int(np.argmax(mi_roll_h))])
roll_top_mi_idx = np.argsort(mi_roll_h)[-min(10, len(mi_roll_h)):][::-1]
roll_top_mi_pairs = [(int(roll_mi_lags[i]), float(mi_roll_h[i])) for i in roll_top_mi_idx]

rng_2b3c = np.random.default_rng(48)
roll_mi_shuffle_runs = 500
shuffle_roll_max_mi = []
for _ in range(roll_mi_shuffle_runs):
    shuffled_h = rng_2b3c.permutation(roll_h)
    shuffled_mi = mutual_information_lags(shuffled_h, roll_mi_max_lag, roll_mi_bins)
    shuffle_roll_max_mi.append(float(shuffled_mi.max()))
shuffle_roll_max_mi = np.asarray(shuffle_roll_max_mi, dtype=float)

shuffle_roll_mi_mean = float(shuffle_roll_max_mi.mean())
shuffle_roll_mi_std = float(shuffle_roll_max_mi.std(ddof=1))
shuffle_roll_mi_p = float(np.mean(shuffle_roll_max_mi >= roll_max_mi))
shuffle_roll_mi_z = (
    (roll_max_mi - shuffle_roll_mi_mean) / (shuffle_roll_mi_std + 1e-12)
)

if shuffle_roll_mi_p <= 0.05:
    roll_mi_note = "rolling H ima MI signal iznad shuffled reference"
else:
    roll_mi_note = "rolling H nema jak MI signal iznad shuffled reference"

print()
print("KORAK 2b3c: Aparat 2b Hurst/R-S + Test 3c Mutual Information")
print(f"  rolling H max MI lag 1..{roll_mi_max_lag}: {roll_max_mi:.6f} bits "
      f"(lag={roll_max_mi_lag})")
print(f"  shuffled max MI: mean={shuffle_roll_mi_mean:.6f} std={shuffle_roll_mi_std:.6f} "
      f"z={shuffle_roll_mi_z:.2f} p={shuffle_roll_mi_p:.4f}")
print(f"  ⇒ {roll_mi_note}")
print()

fig2b3c, ax2b3c = plt.subplots(1, 3, figsize=(16, 5))
fig2b3c.suptitle("KORAK 2b3c: Hurst/R-S aparat + Mutual Information test",
                 fontsize=13, fontweight="bold")

ax2b3c[0].plot(roll_mi_lags, mi_roll_h, "o-", markersize=4, color="darkslateblue")
ax2b3c[0].set_title("MI rolling/local Hurst niza")
ax2b3c[0].set_xlabel("lag")
ax2b3c[0].set_ylabel("MI [bits]")
ax2b3c[0].grid(True, alpha=0.25)

ax2b3c[1].plot(mi_f_control_lags, mi_f_control, "o-", markersize=3, color="steelblue")
ax2b3c[1].set_title("Kontrola: MI f(t)")
ax2b3c[1].set_xlabel("lag")
ax2b3c[1].set_ylabel("MI [bits]")
ax2b3c[1].grid(True, alpha=0.25)

ax2b3c[2].hist(shuffle_roll_max_mi, bins=24, color="lightgray", edgecolor="white")
ax2b3c[2].axvline(roll_max_mi, color="crimson", linewidth=2,
                  label=f"observed={roll_max_mi:.4f}")
ax2b3c[2].axvline(shuffle_roll_mi_mean, color="black", linestyle="--",
                  label=f"shuffle mean={shuffle_roll_mi_mean:.4f}")
ax2b3c[2].set_title("Shuffled rolling H max MI")
ax2b3c[2].set_xlabel("max MI [bits]")
ax2b3c[2].set_ylabel("broj")
ax2b3c[2].legend(fontsize=8)

for a in ax2b3c:
    a.spines["top"].set_visible(False)
    a.spines["right"].set_visible(False)

fig2b3c.tight_layout()
fig2b3c.savefig(PNG_PATH_2B3C, dpi=150, bbox_inches="tight")
plt.show()

with open(TXT_PATH, "a", encoding="utf-8") as f:
    f.write("\n")
    f.write("=" * 60 + "\n")
    f.write("KORAK 2b3c: Aparat 2b Hurst/R-S + Test 3c Mutual Information\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"  PNG:                  {PNG_PATH_2B3C}\n\n")
    f.write("Mutual Information nad rolling/local Hurst nizom:\n")
    f.write(f"  broj rolling H tacaka = {len(roll_h)}\n")
    f.write(f"  max lag               = {roll_mi_max_lag}\n")
    f.write(f"  bins                  = {roll_mi_bins}\n")
    f.write(f"  max MI                = {roll_max_mi:.8f} bits\n")
    f.write(f"  max MI lag            = {roll_max_mi_lag}\n\n")
    f.write("Shuffled rolling H max MI referenca:\n")
    f.write(f"  runs                  = {roll_mi_shuffle_runs}\n")
    f.write(f"  mean                  = {shuffle_roll_mi_mean:.8f}\n")
    f.write(f"  std                   = {shuffle_roll_mi_std:.8f}\n")
    f.write(f"  z                     = {shuffle_roll_mi_z:.8f}\n")
    f.write(f"  p(shuffled >= obs)    = {shuffle_roll_mi_p:.8f}\n")
    f.write(f"  interpret.            = {roll_mi_note}\n\n")
    f.write("Top MI lagovi rolling H:\n")
    f.write(f"  {'lag':<8}{'MI [bits]':>16}\n")
    for lag, val in roll_top_mi_pairs:
        f.write(f"  {lag:<8}{val:>16,.8f}\n")
    f.write("\n")

    elapsed_2b3c = time.time() - T0_2B3C
    f.write(f"Vreme KORAKA 2b3c: {timedelta(seconds=int(elapsed_2b3c))} ({elapsed_2b3c:.1f} s)\n")
    f.write("\nKraj KORAKA 2b3c.\n")

print(f"PNG saved → {PNG_PATH_2B3C}")
print(f"TXT updated → {TXT_PATH}")
print(f"Vreme KORAKA 2b3c: {timedelta(seconds=int(time.time()-T0_2B3C))} "
      f"({time.time()-T0_2B3C:.1f} s)")
print()
"""
KORAK 2b3c (Aparat 2b: Hurst/R-S + Test 3c: Mutual Information)
u 1_KarlWeierstrass_v2.py. Posebna PNG slika, TXT proširen.

KORAK 2b3c — MI nad rolling/local Hurst nizom,
             MI nad f(t) kao kontrola,
             shuffled max MI referenca.
             Crta rezultat u 1_KarlWeierstrass_v2_2b3c.png.


2b3c: Aparat 2b Hurst/R-S + Test 3c Mutual Information. MI nad rolling/local Hurst nizom, MI nad f(t) kao kontrolu i shuffled MI referencu, sa novim PNG _2b3c.

2b3c: MI rolling H režima, MI f(t) kao kontrolu i shuffled max MI referencu za rolling H.

Dodato:
1_KarlWeierstrass_v2_2b3c.png
MI nad rolling/local Hurst nizom
MI nad f(t) kao kontrola
shuffled max MI referenca (500 permutacija)
TXT append i komentar-log blok


KORAK 2b3c: Aparat 2b Hurst/R-S + Test 3c Mutual Information
  rolling H max MI lag 1..12: 1.357025 bits (lag=5)
  shuffled max MI: mean=1.506056 std=0.113304 z=-1.32 p=0.9220
  ⇒ rolling H nema jak MI signal iznad shuffled reference

PNG saved →   /1_KarlWeierstrass_v2_2b3c.png
TXT updated → /1_KarlWeierstrass_v2.txt
Vreme KORAKA 2b3c: 0:00:16 (16.8 s)
"""



# ─────────────────────────────────────────────────────────────────────
# KORAK 2b3d: Aparat 2b Hurst/R-S + Test 3d Entropy
#   Cilj: izmeriti kompleksnost lokalnih Hurst rezima kroz vreme.
#   Rolling H niz je kratak, zato se sample entropy tretira kao indikativna.
# ─────────────────────────────────────────────────────────────────────
T0_2B3D = time.time()

roll_sampen_m = 2
roll_sampen_r = 0.35
roll_sampen, roll_sampen_a, roll_sampen_b, roll_sampen_n = sample_entropy(
    roll_h, m=roll_sampen_m, r=roll_sampen_r, max_points=len(roll_h)
)
f_sampen_control, _, _, f_sampen_n = sample_entropy(
    lex_idx, m=roll_sampen_m, r=0.2, max_points=1200
)

roll_pe_orders = [3, 4, 5]
roll_pe_rows = []
f_pe_rows = []
for order in roll_pe_orders:
    pe, pe_norm, patterns = permutation_entropy(roll_h, order=order, delay=1)
    roll_pe_rows.append((order, pe, pe_norm, patterns, math.factorial(order)))
    pe_f, pe_f_norm, patterns_f = permutation_entropy(lex_idx, order=order, delay=1)
    f_pe_rows.append((order, pe_f, pe_f_norm, patterns_f, math.factorial(order)))

rng_2b3d = np.random.default_rng(49)
roll_entropy_shuffle_runs = 500
shuffle_roll_sampen = []
shuffle_roll_pe4 = []
for _ in range(roll_entropy_shuffle_runs):
    shuffled_h = rng_2b3d.permutation(roll_h)
    se, _, _, _ = sample_entropy(
        shuffled_h, m=roll_sampen_m, r=roll_sampen_r, max_points=len(shuffled_h)
    )
    _, pe4_norm_shuf, _ = permutation_entropy(shuffled_h, order=4, delay=1)
    if np.isfinite(se):
        shuffle_roll_sampen.append(se)
    shuffle_roll_pe4.append(pe4_norm_shuf)

shuffle_roll_sampen = np.asarray(shuffle_roll_sampen, dtype=float)
shuffle_roll_pe4 = np.asarray(shuffle_roll_pe4, dtype=float)
roll_pe4_norm = roll_pe_rows[1][2]

if len(shuffle_roll_sampen) > 0 and np.isfinite(roll_sampen):
    shuffle_roll_sampen_mean = float(shuffle_roll_sampen.mean())
    shuffle_roll_sampen_std = float(shuffle_roll_sampen.std(ddof=1))
    shuffle_roll_sampen_p_low = float(np.mean(shuffle_roll_sampen <= roll_sampen))
    shuffle_roll_sampen_z = (
        (roll_sampen - shuffle_roll_sampen_mean) / (shuffle_roll_sampen_std + 1e-12)
    )
else:
    shuffle_roll_sampen_mean = float("nan")
    shuffle_roll_sampen_std = float("nan")
    shuffle_roll_sampen_p_low = float("nan")
    shuffle_roll_sampen_z = float("nan")

shuffle_roll_pe4_mean = float(shuffle_roll_pe4.mean())
shuffle_roll_pe4_std = float(shuffle_roll_pe4.std(ddof=1))
shuffle_roll_pe4_p_low = float(np.mean(shuffle_roll_pe4 <= roll_pe4_norm))
shuffle_roll_pe4_z = (
    (roll_pe4_norm - shuffle_roll_pe4_mean) / (shuffle_roll_pe4_std + 1e-12)
)

if (np.isfinite(shuffle_roll_sampen_p_low) and shuffle_roll_sampen_p_low <= 0.05) or (
    shuffle_roll_pe4_p_low <= 0.05
):
    roll_entropy_note = "entropy rolling H je niza od shuffled reference (moguca struktura)"
else:
    roll_entropy_note = "entropy rolling H je blizu shuffled reference"

print()
print("KORAK 2b3d: Aparat 2b Hurst/R-S + Test 3d Sample / Permutation entropy")
print(f"  Sample entropy rolling H: {roll_sampen:.4f} "
      f"(m={roll_sampen_m}, r={roll_sampen_r}, n={roll_sampen_n})")
print(f"  Sample entropy f(t) kontrola: {f_sampen_control:.4f} (n={f_sampen_n})")
print(f"  Permutation entropy rolling H order=4: {roll_pe4_norm:.4f} normalizovano")
print(f"  shuffled SampEn: mean={shuffle_roll_sampen_mean:.4f} "
      f"std={shuffle_roll_sampen_std:.4f} z={shuffle_roll_sampen_z:.2f} "
      f"p_low={shuffle_roll_sampen_p_low:.4f}")
print(f"  shuffled PE4: mean={shuffle_roll_pe4_mean:.4f} std={shuffle_roll_pe4_std:.4f} "
      f"z={shuffle_roll_pe4_z:.2f} p_low={shuffle_roll_pe4_p_low:.4f}")
print(f"  ⇒ {roll_entropy_note}")
print()

fig2b3d, ax2b3d = plt.subplots(1, 3, figsize=(16, 5))
fig2b3d.suptitle("KORAK 2b3d: Hurst/R-S aparat + Sample / Permutation entropy",
                 fontsize=13, fontweight="bold")

orders = np.array([row[0] for row in roll_pe_rows], dtype=int)
roll_pe_norms = np.array([row[2] for row in roll_pe_rows], dtype=float)
f_pe_norms = np.array([row[2] for row in f_pe_rows], dtype=float)
ax2b3d[0].plot(orders, roll_pe_norms, "o-", color="darkslateblue", label="rolling H")
ax2b3d[0].plot(orders, f_pe_norms, "o-", color="steelblue", label="f(t)")
ax2b3d[0].set_ylim(0, 1.05)
ax2b3d[0].set_title("Normalizovana permutation entropy")
ax2b3d[0].set_xlabel("order")
ax2b3d[0].set_ylabel("PE / max PE")
ax2b3d[0].legend(fontsize=8)
ax2b3d[0].grid(True, alpha=0.25)

if len(shuffle_roll_sampen) > 0:
    ax2b3d[1].hist(shuffle_roll_sampen, bins=20, color="lightgray", edgecolor="white")
    ax2b3d[1].axvline(roll_sampen, color="crimson", linewidth=2,
                      label=f"observed={roll_sampen:.3f}")
    ax2b3d[1].axvline(shuffle_roll_sampen_mean, color="black", linestyle="--",
                      label=f"shuffle mean={shuffle_roll_sampen_mean:.3f}")
ax2b3d[1].set_title("Shuffled rolling H Sample Entropy")
ax2b3d[1].set_xlabel("Sample entropy")
ax2b3d[1].set_ylabel("broj")
ax2b3d[1].legend(fontsize=8)

ax2b3d[2].hist(shuffle_roll_pe4, bins=20, color="lightgray", edgecolor="white")
ax2b3d[2].axvline(roll_pe4_norm, color="crimson", linewidth=2,
                  label=f"observed={roll_pe4_norm:.3f}")
ax2b3d[2].axvline(shuffle_roll_pe4_mean, color="black", linestyle="--",
                  label=f"shuffle mean={shuffle_roll_pe4_mean:.3f}")
ax2b3d[2].set_title("Shuffled rolling H PE order=4")
ax2b3d[2].set_xlabel("normalizovana PE")
ax2b3d[2].set_ylabel("broj")
ax2b3d[2].legend(fontsize=8)

for a in ax2b3d:
    a.spines["top"].set_visible(False)
    a.spines["right"].set_visible(False)

fig2b3d.tight_layout()
fig2b3d.savefig(PNG_PATH_2B3D, dpi=150, bbox_inches="tight")
plt.show()

with open(TXT_PATH, "a", encoding="utf-8") as f:
    f.write("\n")
    f.write("=" * 60 + "\n")
    f.write("KORAK 2b3d: Aparat 2b Hurst/R-S + Test 3d Sample / Permutation entropy\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"  PNG:                  {PNG_PATH_2B3D}\n\n")
    f.write("Sample entropy rolling/local Hurst niza:\n")
    f.write(f"  m                     = {roll_sampen_m}\n")
    f.write(f"  r                     = {roll_sampen_r}\n")
    f.write(f"  n used                = {roll_sampen_n}\n")
    f.write(f"  SampEn(rolling H)     = {roll_sampen:.8f}\n")
    f.write(f"  SampEn(f(t) kontrola) = {f_sampen_control:.8f}\n")
    f.write(f"  A count               = {roll_sampen_a}\n")
    f.write(f"  B count               = {roll_sampen_b}\n\n")
    f.write("Permutation entropy rolling H:\n")
    f.write(f"  {'order':<8}{'PE bits':>16}{'PE norm':>16}{'patterns':>14}{'max':>10}\n")
    for order, pe, pe_norm, patterns, max_patterns in roll_pe_rows:
        f.write(f"  {order:<8}{pe:>16,.8f}{pe_norm:>16,.8f}{patterns:>14}{max_patterns:>10}\n")
    f.write("\n")
    f.write("Shuffled entropy referenca:\n")
    f.write(f"  runs                  = {roll_entropy_shuffle_runs}\n")
    f.write(f"  SampEn finite runs    = {len(shuffle_roll_sampen)}\n")
    f.write(f"  SampEn mean           = {shuffle_roll_sampen_mean:.8f}\n")
    f.write(f"  SampEn std            = {shuffle_roll_sampen_std:.8f}\n")
    f.write(f"  SampEn z              = {shuffle_roll_sampen_z:.8f}\n")
    f.write(f"  SampEn p_low          = {shuffle_roll_sampen_p_low:.8f}\n")
    f.write(f"  PE4 mean              = {shuffle_roll_pe4_mean:.8f}\n")
    f.write(f"  PE4 std               = {shuffle_roll_pe4_std:.8f}\n")
    f.write(f"  PE4 z                 = {shuffle_roll_pe4_z:.8f}\n")
    f.write(f"  PE4 p_low             = {shuffle_roll_pe4_p_low:.8f}\n")
    f.write(f"  interpret.            = {roll_entropy_note}\n\n")

    elapsed_2b3d = time.time() - T0_2B3D
    f.write(f"Vreme KORAKA 2b3d: {timedelta(seconds=int(elapsed_2b3d))} ({elapsed_2b3d:.1f} s)\n")
    f.write("\nKraj KORAKA 2b3d.\n")

print(f"PNG saved → {PNG_PATH_2B3D}")
print(f"TXT updated → {TXT_PATH}")
print(f"Vreme KORAKA 2b3d: {timedelta(seconds=int(time.time()-T0_2B3D))} "
      f"({time.time()-T0_2B3D:.1f} s)")
print()
"""
KORAK 2b3d (Aparat 2b: Hurst/R-S + Test 3d: Sample / Permutation entropy)
u 1_KarlWeierstrass_v2.py. Posebna PNG slika, TXT proširen.

KORAK 2b3d — sample entropy nad rolling/local Hurst nizom,
             permutation entropy nad rolling/local Hurst nizom,
             shuffled entropy referenca.
             Crta rezultat u 1_KarlWeierstrass_v2_2b3d.png.


2b3d: Aparat 2b Hurst/R-S + Test 3d Sample / Permutation entropy. Radiću entropiju nad rolling/local Hurst nizom, plus f(t) kao kontrolu i shuffled referencu, novi PNG _2b3d.

2b3d: Sample entropy + permutation entropy nad rolling Hurst nizom, f(t) kao kontrola i shuffled entropy referenca za rolling H.

Dodato:
1_KarlWeierstrass_v2_2b3d.png
Sample entropy nad rolling/local Hurst nizom
Permutation entropy nad rolling/local Hurst nizom
f(t) kao kontrola
shuffled entropy referenca (500 permutacija)
TXT append i komentar-log blok


KORAK 2b3d: Aparat 2b Hurst/R-S + Test 3d Sample / Permutation entropy
  Sample entropy rolling H: 2.5649 (m=2, r=0.35, n=31)
  Sample entropy f(t) kontrola: 2.1954 (n=1200)
  Permutation entropy rolling H order=4: 0.8713 normalizovano
  shuffled SampEn: mean=1.8468 std=0.5096 z=1.41 p_low=0.9042
  shuffled PE4: mean=0.8598 std=0.0392 z=0.29 p_low=0.5880
  ⇒ entropy rolling H je blizu shuffled reference

PNG saved →   /1_KarlWeierstrass_v2_2b3d.png
TXT updated → /1_KarlWeierstrass_v2.txt
Vreme KORAKA 2b3d: 0:00:16 (16.7 s)
"""



# ─────────────────────────────────────────────────────────────────────
# KORAK 2b3e: Aparat 2b Hurst/R-S + Test 3e NIST baterija
#   Cilj: proveriti binarnu nasumicnost lokalnih Hurst rezima.
#   Bit = 1 znaci da je lokalni H iznad Brown/random reference H=0.5.
# ─────────────────────────────────────────────────────────────────────
T0_2B3E = time.time()

roll_h_bits = (roll_h > 0.5).astype(int)
roll_h_n = len(roll_h_bits)
roll_h_ones = int(roll_h_bits.sum())
roll_h_zeros = int(roll_h_n - roll_h_ones)

roll_monobit_p, roll_monobit_s = nist_monobit(roll_h_bits)
roll_runs_p, roll_runs_count, roll_runs_pi = nist_runs(roll_h_bits)
roll_block_p, roll_block_chi2, roll_block_count = nist_block_frequency(
    roll_h_bits, block_size=8
)
roll_cusum_p, roll_cusum_z, roll_cusum_walk = nist_cumulative_sums(roll_h_bits)
roll_apen_p, roll_apen_value, roll_apen_chi2, roll_apen_df = nist_approximate_entropy(
    roll_h_bits, m=2
)

roll_nist_rows = [
    ("Monobit frequency", roll_monobit_p, roll_monobit_s),
    ("Runs", roll_runs_p, roll_runs_count),
    ("Block frequency", roll_block_p, roll_block_chi2),
    ("Cumulative sums", roll_cusum_p, roll_cusum_z),
    ("Approx entropy", roll_apen_p, roll_apen_value),
]
roll_nist_pass_count = int(sum(p > 0.05 for _, p, _ in roll_nist_rows if np.isfinite(p)))
roll_nist_total = int(sum(np.isfinite(p) for _, p, _ in roll_nist_rows))

f_control_bits = nist_bits_from_series(lex_idx)
f_monobit_p, f_monobit_s = nist_monobit(f_control_bits)
f_runs_p, f_runs_count, f_runs_pi = nist_runs(f_control_bits)
f_block_p, f_block_chi2, f_block_count = nist_block_frequency(f_control_bits, block_size=128)
f_cusum_p, f_cusum_z, _ = nist_cumulative_sums(f_control_bits)
f_apen_p, f_apen_value, f_apen_chi2, f_apen_df = nist_approximate_entropy(f_control_bits, m=2)
f_nist_rows = [
    ("Monobit frequency", f_monobit_p, f_monobit_s),
    ("Runs", f_runs_p, f_runs_count),
    ("Block frequency", f_block_p, f_block_chi2),
    ("Cumulative sums", f_cusum_p, f_cusum_z),
    ("Approx entropy", f_apen_p, f_apen_value),
]
f_nist_pass_count = int(sum(p > 0.05 for _, p, _ in f_nist_rows if np.isfinite(p)))
f_nist_total = int(sum(np.isfinite(p) for _, p, _ in f_nist_rows))

if roll_nist_pass_count == roll_nist_total:
    roll_nist_note = "rolling H bitovi prolaze sve NIST-style testove"
elif roll_nist_pass_count >= max(1, roll_nist_total - 1):
    roll_nist_note = "rolling H bitovi uglavnom prolaze, slab signal za proveru"
else:
    roll_nist_note = "rolling H bitovi padaju vise NIST-style testova"

print()
print("KORAK 2b3e: Aparat 2b Hurst/R-S + Test 3e NIST baterija")
print(f"  rolling H bits: n={roll_h_n}  zeros={roll_h_zeros}  ones={roll_h_ones}")
for name, p, stat_val in roll_nist_rows:
    print(f"  {name:<20} p={p:.6f}  stat={stat_val}")
print(f"  prolaz rolling H: {roll_nist_pass_count}/{roll_nist_total}  ⇒ {roll_nist_note}")
print(f"  kontrola f(t) prolaz: {f_nist_pass_count}/{f_nist_total}")
print()

fig2b3e, ax2b3e = plt.subplots(1, 3, figsize=(16, 5))
fig2b3e.suptitle("KORAK 2b3e: Hurst/R-S aparat + NIST-style testovi",
                 fontsize=13, fontweight="bold")

ax2b3e[0].bar(["H<=0.5", "H>0.5"], [roll_h_zeros, roll_h_ones],
              color=["steelblue", "darkslateblue"])
ax2b3e[0].set_title("Binarizovani rolling H rezimi")
ax2b3e[0].set_xlabel("bit")
ax2b3e[0].set_ylabel("broj")
ax2b3e[0].grid(True, alpha=0.2, axis="y")

names = [row[0] for row in roll_nist_rows]
roll_pvals = np.array([row[1] for row in roll_nist_rows], dtype=float)
colors = ["seagreen" if p > 0.05 else "crimson" for p in roll_pvals]
ax2b3e[1].barh(names, roll_pvals, color=colors)
ax2b3e[1].axvline(0.05, color="black", linestyle="--", linewidth=1.2)
ax2b3e[1].set_xlim(0, 1)
ax2b3e[1].set_title("Rolling H NIST-style p-vrednosti")
ax2b3e[1].set_xlabel("p-value")
ax2b3e[1].grid(True, alpha=0.2, axis="x")

ax2b3e[2].plot(np.arange(1, roll_h_n + 1), roll_cusum_walk,
               linewidth=1.2, marker="o", markersize=3, color="purple")
ax2b3e[2].axhline(0, color="black", linewidth=0.6)
ax2b3e[2].set_title(f"Rolling H cumulative sums (z={roll_cusum_z})")
ax2b3e[2].set_xlabel("rolling window index")
ax2b3e[2].set_ylabel("cum sum")
ax2b3e[2].grid(True, alpha=0.25)

for a in ax2b3e:
    a.spines["top"].set_visible(False)
    a.spines["right"].set_visible(False)

fig2b3e.tight_layout()
fig2b3e.savefig(PNG_PATH_2B3E, dpi=150, bbox_inches="tight")
plt.show()

with open(TXT_PATH, "a", encoding="utf-8") as f:
    f.write("\n")
    f.write("=" * 60 + "\n")
    f.write("KORAK 2b3e: Aparat 2b Hurst/R-S + Test 3e NIST baterija\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"  PNG:                  {PNG_PATH_2B3E}\n\n")
    f.write("Binarizacija rolling/local Hurst niza:\n")
    f.write("  bit = 1 ako je lokalni H > 0.5, inace 0\n")
    f.write(f"  n bits                = {roll_h_n}\n")
    f.write(f"  zeros                 = {roll_h_zeros}\n")
    f.write(f"  ones                  = {roll_h_ones}\n\n")
    f.write("NIST-style testovi nad rolling H bitovima (prolaz ako p > 0.05):\n")
    f.write(f"  {'test':<22}{'p-value':>14}{'stat':>18}{'pass':>10}\n")
    for name, p, stat_val in roll_nist_rows:
        f.write(f"  {name:<22}{p:>14,.8f}{float(stat_val):>18,.8f}{str(p > 0.05):>10}\n")
    f.write("\n")
    f.write("Detalji rolling H:\n")
    f.write(f"  Runs pi               = {roll_runs_pi:.8f}\n")
    f.write(f"  Block size            = 8\n")
    f.write(f"  Block count           = {roll_block_count}\n")
    f.write(f"  Approx entropy m      = 2\n")
    f.write(f"  Approx entropy chi2   = {roll_apen_chi2:.8f}\n")
    f.write(f"  Approx entropy df     = {roll_apen_df}\n")
    f.write(f"  pass count            = {roll_nist_pass_count}/{roll_nist_total}\n")
    f.write(f"  interpret.            = {roll_nist_note}\n\n")
    f.write("Kontrola: NIST-style testovi nad f(t) binarizacijom:\n")
    f.write(f"  {'test':<22}{'p-value':>14}{'stat':>18}{'pass':>10}\n")
    for name, p, stat_val in f_nist_rows:
        f.write(f"  {name:<22}{p:>14,.8f}{float(stat_val):>18,.8f}{str(p > 0.05):>10}\n")
    f.write(f"  pass count            = {f_nist_pass_count}/{f_nist_total}\n")
    f.write(f"  f(t) runs pi          = {f_runs_pi:.8f}\n")
    f.write(f"  f(t) block count      = {f_block_count}\n")
    f.write(f"  f(t) approx chi2      = {f_apen_chi2:.8f}\n")
    f.write(f"  f(t) approx df        = {f_apen_df}\n\n")

    elapsed_2b3e = time.time() - T0_2B3E
    f.write(f"Vreme KORAKA 2b3e: {timedelta(seconds=int(elapsed_2b3e))} ({elapsed_2b3e:.1f} s)\n")
    f.write("\nKraj KORAKA 2b3e.\n")

print(f"PNG saved → {PNG_PATH_2B3E}")
print(f"TXT updated → {TXT_PATH}")
print(f"Vreme KORAKA 2b3e: {timedelta(seconds=int(time.time()-T0_2B3E))} "
      f"({time.time()-T0_2B3E:.1f} s)")
print()
"""
KORAK 2b3e (Aparat 2b: Hurst/R-S + Test 3e: NIST baterija)
u 1_KarlWeierstrass_v2.py. Posebna PNG slika, TXT proširen.

KORAK 2b3e — NIST-style testovi nad rolling H > 0.5 bitovima:
             monobit, runs, block frequency, cumulative sums,
             approximate entropy; f(t) binarizacija kao kontrola.
             Crta rezultat u 1_KarlWeierstrass_v2_2b3e.png.


2b3e: Aparat 2b Hurst/R-S + Test 3e NIST baterija. Radiću NIST-style testove nad binarizovanim rolling/local Hurst nizom, uz f(t) kao kontrolu, novi PNG _2b3e, TXT append.

2b3e: NIST-style testovi nad bitovima rolling H > 0.5, plus kontrolni NIST nad f(t) binarizacijom.

Dodato:
1_KarlWeierstrass_v2_2b3e.png
NIST-style testovi nad rolling H > 0.5 bitovima
f(t) NIST binarizacija kao kontrola
Monobit, Runs, Block frequency, Cumulative sums, Approximate entropy
TXT append i komentar-log blok


KORAK 2b3e: Aparat 2b Hurst/R-S + Test 3e NIST baterija
  rolling H bits: n=31  zeros=0  ones=31
  Monobit frequency    p=0.000000  stat=5.567764362830022
  Runs                 p=0.000000  stat=0
  Block frequency      p=0.000025  stat=24.0
  Cumulative sums      p=0.000000  stat=31
  Approx entropy       p=0.000000  stat=0.0
  prolaz rolling H: 0/5  ⇒ rolling H bitovi padaju vise NIST-style testova
  kontrola f(t) prolaz: 4/5

PNG saved →   /1_KarlWeierstrass_v2_2b3e.png
TXT updated → /1_KarlWeierstrass_v2.txt
Vreme KORAKA 2b3e: 0:00:15 (15.8 s)
"""



# ─────────────────────────────────────────────────────────────────────
# KORAK 2c3a: Aparat 2c Fraktalna dimenzija + Test 3a Hurst eksponent
#   Cilj: proveriti vezu izmedju fraktalne dimenzije i Hurst eksponenta.
#   Za samoslicne grafove cesto vazi gruba veza D ~= 2 - H.
# ─────────────────────────────────────────────────────────────────────
T0_2C3A = time.time()

fd_hurst_reference = 2.0 - hurst_f
fd_hurst_gap = higuchi_f - fd_hurst_reference

fd_roll_window = 768
fd_roll_step = 128
fd_roll_centers, fd_roll, fd_roll_r2 = rolling_higuchi_fd(
    lex_idx, window=fd_roll_window, step=fd_roll_step, kmax=32
)
fd_from_roll_h = 2.0 - roll_h
fd_roll_gap = fd_roll - fd_from_roll_h

rng_2c3a = np.random.default_rng(50)
fd_shuffle_runs = 100
shuffle_fd = []
for _ in range(fd_shuffle_runs):
    shuffled_f = rng_2c3a.permutation(lex_idx)
    fd_shuf, _, _, _, _ = higuchi_fd(normalize01(shuffled_f), kmax=64)
    shuffle_fd.append(fd_shuf)
shuffle_fd = np.asarray(shuffle_fd, dtype=float)

shuffle_fd_mean = float(shuffle_fd.mean())
shuffle_fd_std = float(shuffle_fd.std(ddof=1))
shuffle_fd_p_high = float(np.mean(shuffle_fd >= higuchi_f))
shuffle_fd_p_low = float(np.mean(shuffle_fd <= higuchi_f))
shuffle_fd_z = (higuchi_f - shuffle_fd_mean) / (shuffle_fd_std + 1e-12)

fd_roll_mean = float(fd_roll.mean())
fd_roll_std = float(fd_roll.std(ddof=1))
fd_roll_gap_mean = float(fd_roll_gap.mean())
fd_roll_gap_std = float(fd_roll_gap.std(ddof=1))
corr_fd_h = float(np.corrcoef(fd_roll, roll_h)[0, 1])

if shuffle_fd_p_high <= 0.05 and higuchi_f > shuffle_fd_mean:
    fd_hurst_note = "globalni FD je iznad shuffled reference"
elif shuffle_fd_p_low <= 0.05 and higuchi_f < shuffle_fd_mean:
    fd_hurst_note = "globalni FD je ispod shuffled reference"
else:
    fd_hurst_note = "globalni FD nije jak odmak od shuffled reference"

print()
print("KORAK 2c3a: Aparat 2c Fraktalna dimenzija + Test 3a Hurst eksponent")
print(f"  Higuchi FD(f(t)) = {higuchi_f:.4f}   R²={higuchi_r2:.4f}")
print(f"  Hurst referenca 2-H = {fd_hurst_reference:.4f}   gap FD-(2-H)={fd_hurst_gap:.4f}")
print(f"  rolling FD: mean={fd_roll_mean:.4f} std={fd_roll_std:.4f}")
print(f"  rolling gap FD-(2-H): mean={fd_roll_gap_mean:.4f} std={fd_roll_gap_std:.4f}")
print(f"  corr(rolling FD, rolling H) = {corr_fd_h:.4f}")
print(f"  shuffled FD: mean={shuffle_fd_mean:.4f} std={shuffle_fd_std:.4f} "
      f"z={shuffle_fd_z:.2f} p_high={shuffle_fd_p_high:.4f}")
print(f"  ⇒ {fd_hurst_note}")
print()

fig2c3a, ax2c3a = plt.subplots(1, 3, figsize=(16, 5))
fig2c3a.suptitle("KORAK 2c3a: Fraktalna dimenzija + Hurst test",
                 fontsize=13, fontweight="bold")

ax2c3a[0].plot(fd_roll_centers, fd_roll, "o-", markersize=3, color="purple",
               label="rolling Higuchi FD")
ax2c3a[0].plot(roll_centers, fd_from_roll_h, "o-", markersize=3, color="darkslateblue",
               label="2 - rolling H")
ax2c3a[0].set_title("Rolling FD vs Hurst referenca")
ax2c3a[0].set_xlabel("t centar prozora")
ax2c3a[0].set_ylabel("vrednost")
ax2c3a[0].legend(fontsize=8)
ax2c3a[0].grid(True, alpha=0.25)

ax2c3a[1].scatter(roll_h, fd_roll, color="purple", alpha=0.8)
xs_h = np.linspace(float(roll_h.min()), float(roll_h.max()), 100)
ax2c3a[1].plot(xs_h, 2.0 - xs_h, "k--", label="D = 2 - H")
ax2c3a[1].set_title(f"Rolling FD prema rolling H (corr={corr_fd_h:.3f})")
ax2c3a[1].set_xlabel("rolling H")
ax2c3a[1].set_ylabel("rolling Higuchi FD")
ax2c3a[1].legend(fontsize=8)
ax2c3a[1].grid(True, alpha=0.25)

ax2c3a[2].hist(shuffle_fd, bins=22, color="lightgray", edgecolor="white")
ax2c3a[2].axvline(higuchi_f, color="crimson", linewidth=2,
                  label=f"observed FD={higuchi_f:.3f}")
ax2c3a[2].axvline(shuffle_fd_mean, color="black", linestyle="--",
                  label=f"shuffle mean={shuffle_fd_mean:.3f}")
ax2c3a[2].set_title("Shuffled Higuchi FD referenca")
ax2c3a[2].set_xlabel("Higuchi FD shuffled f(t)")
ax2c3a[2].set_ylabel("broj")
ax2c3a[2].legend(fontsize=8)

for a in ax2c3a:
    a.spines["top"].set_visible(False)
    a.spines["right"].set_visible(False)

fig2c3a.tight_layout()
fig2c3a.savefig(PNG_PATH_2C3A, dpi=150, bbox_inches="tight")
plt.show()

with open(TXT_PATH, "a", encoding="utf-8") as f:
    f.write("\n")
    f.write("=" * 60 + "\n")
    f.write("KORAK 2c3a: Aparat 2c Fraktalna dimenzija + Test 3a Hurst eksponent\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"  PNG:                  {PNG_PATH_2C3A}\n\n")
    f.write("Globalna veza FD i Hurst-a:\n")
    f.write(f"  Higuchi FD(f(t))      = {higuchi_f:.8f}\n")
    f.write(f"  Higuchi R^2           = {higuchi_r2:.8f}\n")
    f.write(f"  H(f(t))               = {hurst_f:.8f}\n")
    f.write(f"  2 - H                 = {fd_hurst_reference:.8f}\n")
    f.write(f"  FD - (2-H)            = {fd_hurst_gap:.8f}\n\n")
    f.write("Rolling/local FD prema rolling H:\n")
    f.write(f"  window                = {fd_roll_window}\n")
    f.write(f"  step                  = {fd_roll_step}\n")
    f.write(f"  broj prozora          = {len(fd_roll)}\n")
    f.write(f"  mean FD               = {fd_roll_mean:.8f}\n")
    f.write(f"  std FD                = {fd_roll_std:.8f}\n")
    f.write(f"  mean gap FD-(2-H)     = {fd_roll_gap_mean:.8f}\n")
    f.write(f"  std gap FD-(2-H)      = {fd_roll_gap_std:.8f}\n")
    f.write(f"  corr(FD,H)            = {corr_fd_h:.8f}\n\n")
    f.write("Shuffled FD referenca:\n")
    f.write(f"  runs                  = {fd_shuffle_runs}\n")
    f.write(f"  mean                  = {shuffle_fd_mean:.8f}\n")
    f.write(f"  std                   = {shuffle_fd_std:.8f}\n")
    f.write(f"  z                     = {shuffle_fd_z:.8f}\n")
    f.write(f"  p_high                = {shuffle_fd_p_high:.8f}\n")
    f.write(f"  p_low                 = {shuffle_fd_p_low:.8f}\n")
    f.write(f"  interpret.            = {fd_hurst_note}\n\n")
    f.write("Rolling FD tacke:\n")
    f.write(f"  {'center':<10}{'FD':>16}{'R^2':>16}{'2-H':>16}{'gap':>16}\n")
    for center, fdv, r2v, refv, gapv in zip(
        fd_roll_centers.astype(int), fd_roll, fd_roll_r2, fd_from_roll_h, fd_roll_gap
    ):
        f.write(f"  {center:<10}{fdv:>16,.8f}{r2v:>16,.8f}{refv:>16,.8f}{gapv:>16,.8f}\n")
    f.write("\n")

    elapsed_2c3a = time.time() - T0_2C3A
    f.write(f"Vreme KORAKA 2c3a: {timedelta(seconds=int(elapsed_2c3a))} ({elapsed_2c3a:.1f} s)\n")
    f.write("\nKraj KORAKA 2c3a.\n")

print(f"PNG saved → {PNG_PATH_2C3A}")
print(f"TXT updated → {TXT_PATH}")
print(f"Vreme KORAKA 2c3a: {timedelta(seconds=int(time.time()-T0_2C3A))} "
      f"({time.time()-T0_2C3A:.1f} s)")
print()
print("KRAJ KORAKA 2c3a.")
print()
"""
KORAK 2c3a (Aparat 2c: Fraktalna dimenzija + Test 3a: Hurst eksponent)
u 1_KarlWeierstrass_v2.py. Posebna PNG slika, TXT proširen.

KORAK 2c3a — veza Higuchi FD i Hurst reference D ~= 2 - H,
             rolling/local FD kroz vreme,
             shuffled FD referenca.
             Crta rezultat u 1_KarlWeierstrass_v2_2c3a.png.


2c3a: Aparat 2c Fraktalna dimenzija + Test 3a Hurst eksponent. Radiću vezu između Higuchi FD i Hurst reference (D ≈ 2 - H), rolling/local FD kroz vreme i shuffled referencu, sa novim PNG _2c3a.

Dodato:
1_KarlWeierstrass_v2_2c3a.png
veza Higuchi FD i Hurst reference D ~= 2 - H
rolling/local Higuchi FD kroz vreme
poređenje rolling FD sa rolling H
shuffled FD referenca (100 permutacija)
TXT append i komentar-log blok


KORAK 2c3a: Aparat 2c Fraktalna dimenzija + Test 3a Hurst eksponent
  Higuchi FD(f(t)) = 1.9988   R²=1.0000
  Hurst referenca 2-H = 1.4069   gap FD-(2-H)=0.5919
  rolling FD: mean=1.9975 std=0.0028
  rolling gap FD-(2-H): mean=0.5950 std=0.0242
  corr(rolling FD, rolling H) = -0.4357
  shuffled FD: mean=2.0001 std=0.0011 z=-1.11 p_high=0.8500
  ⇒ globalni FD nije jak odmak od shuffled reference

PNG saved →   /1_KarlWeierstrass_v2_2c3a.png
TXT updated → /1_KarlWeierstrass_v2.txt
Vreme KORAKA 2c3a: 0:00:15 (15.5 s)
"""



# ─────────────────────────────────────────────────────────────────────
# KORAK 2c3b: Aparat 2c Fraktalna dimenzija + Test 3b Autokorelacija (ACF)
#   Cilj: proveriti da li lokalna fraktalna dimenzija ima linearno pamcenje.
#   ACF rolling FD niza meri da li se hrapavost grupise kroz vreme.
# ─────────────────────────────────────────────────────────────────────
T0_2C3B = time.time()

fd_acf_max_lag = min(20, max(1, len(fd_roll) - 2))
fd_acf_lags = np.arange(0, fd_acf_max_lag + 1)
acf_fd_roll = autocorr_values(fd_roll, fd_acf_max_lag)
acf_fd_control = autocorr_values(lex_idx, 60)

fd_acf_band = 1.96 / np.sqrt(len(fd_roll))
fd_acf_body = acf_fd_roll[1:]
fd_max_abs_acf = float(np.max(np.abs(fd_acf_body)))
fd_sig_lag_count = int(np.sum(np.abs(fd_acf_body) > fd_acf_band))
fd_top_idx = np.argsort(np.abs(fd_acf_body))[-min(10, len(fd_acf_body)):][::-1] + 1
fd_top_acf_pairs = [(int(lag), float(acf_fd_roll[lag])) for lag in fd_top_idx]

fd_lb_h = min(10, fd_acf_max_lag)
fd_lb_q, fd_lb_p = ljung_box_approx(acf_fd_roll, len(fd_roll), fd_lb_h)

rng_2c3b = np.random.default_rng(51)
fd_acf_shuffle_runs = 500
shuffle_fd_max_abs_acf = []
for _ in range(fd_acf_shuffle_runs):
    shuffled_fd = rng_2c3b.permutation(fd_roll)
    shuffled_acf = autocorr_values(shuffled_fd, fd_acf_max_lag)
    shuffle_fd_max_abs_acf.append(float(np.max(np.abs(shuffled_acf[1:]))))
shuffle_fd_max_abs_acf = np.asarray(shuffle_fd_max_abs_acf, dtype=float)
shuffle_fd_acf_mean = float(shuffle_fd_max_abs_acf.mean())
shuffle_fd_acf_std = float(shuffle_fd_max_abs_acf.std(ddof=1))
shuffle_fd_acf_p = float(np.mean(shuffle_fd_max_abs_acf >= fd_max_abs_acf))
shuffle_fd_acf_z = (fd_max_abs_acf - shuffle_fd_acf_mean) / (shuffle_fd_acf_std + 1e-12)

if fd_lb_p <= 0.05 or shuffle_fd_acf_p <= 0.05:
    fd_acf_note = "rolling FD ima ACF signal iznad shuffled reference"
else:
    fd_acf_note = "rolling FD nema jak ACF signal iznad shuffled reference"

print()
print("KORAK 2c3b: Aparat 2c Fraktalna dimenzija + Test 3b Autokorelacija (ACF)")
print(f"  rolling FD max |ACF| lag 1..{fd_acf_max_lag}: {fd_max_abs_acf:.4f}")
print(f"  95% band: +/-{fd_acf_band:.4f}   znacajnih lagova: "
      f"{fd_sig_lag_count}/{fd_acf_max_lag}")
print(f"  Ljung-Box aproks. h={fd_lb_h}: Q={fd_lb_q:.2f}  p={fd_lb_p:.4f}")
print(f"  shuffled max|ACF|: mean={shuffle_fd_acf_mean:.4f} std={shuffle_fd_acf_std:.4f} "
      f"z={shuffle_fd_acf_z:.2f} p={shuffle_fd_acf_p:.4f}")
print(f"  ⇒ {fd_acf_note}")
print()

fig2c3b, ax2c3b = plt.subplots(1, 3, figsize=(16, 5))
fig2c3b.suptitle("KORAK 2c3b: Fraktalna dimenzija + ACF test",
                 fontsize=13, fontweight="bold")

ax2c3b[0].bar(fd_acf_lags[1:], acf_fd_roll[1:], width=0.8, color="purple")
ax2c3b[0].axhline(fd_acf_band, color="crimson", linestyle="--", linewidth=1.2)
ax2c3b[0].axhline(-fd_acf_band, color="crimson", linestyle="--", linewidth=1.2)
ax2c3b[0].axhline(0, color="black", linewidth=0.6)
ax2c3b[0].set_title("ACF rolling/local FD niza")
ax2c3b[0].set_xlabel("lag")
ax2c3b[0].set_ylabel("ACF")

ax2c3b[1].bar(np.arange(1, 61), acf_fd_control[1:], width=0.8, color="steelblue")
ax2c3b[1].axhline(1.96 / np.sqrt(N), color="crimson", linestyle="--", linewidth=1.2)
ax2c3b[1].axhline(-1.96 / np.sqrt(N), color="crimson", linestyle="--", linewidth=1.2)
ax2c3b[1].axhline(0, color="black", linewidth=0.6)
ax2c3b[1].set_title("Kontrola: ACF f(t)")
ax2c3b[1].set_xlabel("lag")
ax2c3b[1].set_ylabel("ACF")

ax2c3b[2].hist(shuffle_fd_max_abs_acf, bins=24, color="lightgray", edgecolor="white")
ax2c3b[2].axvline(fd_max_abs_acf, color="crimson", linewidth=2,
                  label=f"observed={fd_max_abs_acf:.3f}")
ax2c3b[2].axvline(shuffle_fd_acf_mean, color="black", linestyle="--",
                  label=f"shuffle mean={shuffle_fd_acf_mean:.3f}")
ax2c3b[2].set_title("Shuffled rolling FD max |ACF|")
ax2c3b[2].set_xlabel("max |ACF|")
ax2c3b[2].set_ylabel("broj")
ax2c3b[2].legend(fontsize=8)

for a in ax2c3b:
    a.spines["top"].set_visible(False)
    a.spines["right"].set_visible(False)
    a.grid(True, alpha=0.2)

fig2c3b.tight_layout()
fig2c3b.savefig(PNG_PATH_2C3B, dpi=150, bbox_inches="tight")
plt.show()

with open(TXT_PATH, "a", encoding="utf-8") as f:
    f.write("\n")
    f.write("=" * 60 + "\n")
    f.write("KORAK 2c3b: Aparat 2c Fraktalna dimenzija + Test 3b Autokorelacija (ACF)\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"  PNG:                  {PNG_PATH_2C3B}\n\n")
    f.write("ACF nad rolling/local FD nizom:\n")
    f.write(f"  broj rolling FD tacaka = {len(fd_roll)}\n")
    f.write(f"  max lag               = {fd_acf_max_lag}\n")
    f.write(f"  95% band              = +/-{fd_acf_band:.6f}\n")
    f.write(f"  max |ACF|             = {fd_max_abs_acf:.6f}\n")
    f.write(f"  znacajnih lagova      = {fd_sig_lag_count}/{fd_acf_max_lag}\n")
    f.write(f"  Ljung-Box h           = {fd_lb_h}\n")
    f.write(f"  Ljung-Box Q           = {fd_lb_q:.6f}\n")
    f.write(f"  Ljung-Box p           = {fd_lb_p:.6f}\n\n")
    f.write("Shuffled rolling FD max |ACF| referenca:\n")
    f.write(f"  runs                  = {fd_acf_shuffle_runs}\n")
    f.write(f"  mean                  = {shuffle_fd_acf_mean:.6f}\n")
    f.write(f"  std                   = {shuffle_fd_acf_std:.6f}\n")
    f.write(f"  z                     = {shuffle_fd_acf_z:.6f}\n")
    f.write(f"  p(shuffled >= obs)    = {shuffle_fd_acf_p:.6f}\n")
    f.write(f"  interpret.            = {fd_acf_note}\n\n")
    f.write("Top ACF lagovi rolling FD po apsolutnoj vrednosti:\n")
    f.write(f"  {'lag':<8}{'ACF':>16}\n")
    for lag, val in fd_top_acf_pairs:
        f.write(f"  {lag:<8}{val:>16,.8f}\n")
    f.write("\n")

    elapsed_2c3b = time.time() - T0_2C3B
    f.write(f"Vreme KORAKA 2c3b: {timedelta(seconds=int(elapsed_2c3b))} ({elapsed_2c3b:.1f} s)\n")
    f.write("\nKraj KORAKA 2c3b.\n")

print(f"PNG saved → {PNG_PATH_2C3B}")
print(f"TXT updated → {TXT_PATH}")
print(f"Vreme KORAKA 2c3b: {timedelta(seconds=int(time.time()-T0_2C3B))} "
      f"({time.time()-T0_2C3B:.1f} s)")
print()
print("KRAJ KORAKA 2c3b.")
print()
"""
KORAK 2c3b (Aparat 2c: Fraktalna dimenzija + Test 3b: Autokorelacija / ACF)
u 1_KarlWeierstrass_v2.py. Posebna PNG slika, TXT proširen.

KORAK 2c3b — ACF nad rolling/local FD nizom,
             ACF nad f(t) kao kontrola,
             Ljung-Box i shuffled max|ACF| referenca.
             Crta rezultat u 1_KarlWeierstrass_v2_2c3b.png.


2c3b: Aparat 2c Fraktalna dimenzija + Test 3b Autokorelacija (ACF). Radiću ACF nad rolling/local FD nizom, f(t) kao kontrolu i shuffled ACF referencu, sa novim PNG _2c3b.

2c3b: ACF rolling FD režima, ACF f(t) kao kontrolu, Ljung-Box i shuffled max|ACF| referencu za rolling FD.

Dodato:
1_KarlWeierstrass_v2_2c3b.png
ACF nad rolling/local FD nizom
ACF nad f(t) kao kontrola
Ljung-Box aproksimacija
shuffled max |ACF| referenca (500 permutacija)
TXT append i komentar-log blok


KORAK 2c3b: Aparat 2c Fraktalna dimenzija + Test 3b Autokorelacija (ACF)
  rolling FD max |ACF| lag 1..20: 0.5257
  95% band: +/-0.3520   znacajnih lagova: 2/20
  Ljung-Box aproks. h=10: Q=30.45  p=0.0007
  shuffled max|ACF|: mean=0.3119 std=0.0740 z=2.89 p=0.0140
  ⇒ rolling FD ima ACF signal iznad shuffled reference

PNG saved →   /1_KarlWeierstrass_v2_2c3b.png
TXT updated → /1_KarlWeierstrass_v2.txt
Vreme KORAKA 2c3b: 0:00:14 (14.3 s)
"""



# ─────────────────────────────────────────────────────────────────────
# KORAK 2c3c: Aparat 2c Fraktalna dimenzija + Test 3c Mutual Information
#   Cilj: proveriti nelinearnu zavisnost lokalne fraktalne dimenzije.
#   MI nad rolling FD hvata veze koje ACF ne mora da vidi.
# ─────────────────────────────────────────────────────────────────────
T0_2C3C = time.time()

fd_mi_max_lag = min(12, max(1, len(fd_roll) // 2))
fd_mi_bins = min(8, max(3, len(fd_roll) // 4))
fd_mi_lags = np.arange(1, fd_mi_max_lag + 1)
mi_fd_roll = mutual_information_lags(fd_roll, fd_mi_max_lag, fd_mi_bins)

mi_fd_control_max_lag = 60
mi_fd_control_bins = 16
mi_fd_control = mutual_information_lags(lex_idx, mi_fd_control_max_lag, mi_fd_control_bins)
mi_fd_control_lags = np.arange(1, mi_fd_control_max_lag + 1)

fd_max_mi = float(mi_fd_roll.max())
fd_max_mi_lag = int(fd_mi_lags[int(np.argmax(mi_fd_roll))])
fd_top_mi_idx = np.argsort(mi_fd_roll)[-min(10, len(mi_fd_roll)):][::-1]
fd_top_mi_pairs = [(int(fd_mi_lags[i]), float(mi_fd_roll[i])) for i in fd_top_mi_idx]

rng_2c3c = np.random.default_rng(52)
fd_mi_shuffle_runs = 500
shuffle_fd_max_mi = []
for _ in range(fd_mi_shuffle_runs):
    shuffled_fd = rng_2c3c.permutation(fd_roll)
    shuffled_mi = mutual_information_lags(shuffled_fd, fd_mi_max_lag, fd_mi_bins)
    shuffle_fd_max_mi.append(float(shuffled_mi.max()))
shuffle_fd_max_mi = np.asarray(shuffle_fd_max_mi, dtype=float)

shuffle_fd_mi_mean = float(shuffle_fd_max_mi.mean())
shuffle_fd_mi_std = float(shuffle_fd_max_mi.std(ddof=1))
shuffle_fd_mi_p = float(np.mean(shuffle_fd_max_mi >= fd_max_mi))
shuffle_fd_mi_z = (fd_max_mi - shuffle_fd_mi_mean) / (shuffle_fd_mi_std + 1e-12)

if shuffle_fd_mi_p <= 0.05:
    fd_mi_note = "rolling FD ima MI signal iznad shuffled reference"
else:
    fd_mi_note = "rolling FD nema jak MI signal iznad shuffled reference"

print()
print("KORAK 2c3c: Aparat 2c Fraktalna dimenzija + Test 3c Mutual Information")
print(f"  rolling FD max MI lag 1..{fd_mi_max_lag}: {fd_max_mi:.6f} bits "
      f"(lag={fd_max_mi_lag})")
print(f"  shuffled max MI: mean={shuffle_fd_mi_mean:.6f} std={shuffle_fd_mi_std:.6f} "
      f"z={shuffle_fd_mi_z:.2f} p={shuffle_fd_mi_p:.4f}")
print(f"  ⇒ {fd_mi_note}")
print()

fig2c3c, ax2c3c = plt.subplots(1, 3, figsize=(16, 5))
fig2c3c.suptitle("KORAK 2c3c: Fraktalna dimenzija + Mutual Information test",
                 fontsize=13, fontweight="bold")

ax2c3c[0].plot(fd_mi_lags, mi_fd_roll, "o-", markersize=4, color="purple")
ax2c3c[0].set_title("MI rolling/local FD niza")
ax2c3c[0].set_xlabel("lag")
ax2c3c[0].set_ylabel("MI [bits]")
ax2c3c[0].grid(True, alpha=0.25)

ax2c3c[1].plot(mi_fd_control_lags, mi_fd_control, "o-", markersize=3, color="steelblue")
ax2c3c[1].set_title("Kontrola: MI f(t)")
ax2c3c[1].set_xlabel("lag")
ax2c3c[1].set_ylabel("MI [bits]")
ax2c3c[1].grid(True, alpha=0.25)

ax2c3c[2].hist(shuffle_fd_max_mi, bins=24, color="lightgray", edgecolor="white")
ax2c3c[2].axvline(fd_max_mi, color="crimson", linewidth=2,
                  label=f"observed={fd_max_mi:.4f}")
ax2c3c[2].axvline(shuffle_fd_mi_mean, color="black", linestyle="--",
                  label=f"shuffle mean={shuffle_fd_mi_mean:.4f}")
ax2c3c[2].set_title("Shuffled rolling FD max MI")
ax2c3c[2].set_xlabel("max MI [bits]")
ax2c3c[2].set_ylabel("broj")
ax2c3c[2].legend(fontsize=8)

for a in ax2c3c:
    a.spines["top"].set_visible(False)
    a.spines["right"].set_visible(False)

fig2c3c.tight_layout()
fig2c3c.savefig(PNG_PATH_2C3C, dpi=150, bbox_inches="tight")
plt.show()

with open(TXT_PATH, "a", encoding="utf-8") as f:
    f.write("\n")
    f.write("=" * 60 + "\n")
    f.write("KORAK 2c3c: Aparat 2c Fraktalna dimenzija + Test 3c Mutual Information\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"  PNG:                  {PNG_PATH_2C3C}\n\n")
    f.write("Mutual Information nad rolling/local FD nizom:\n")
    f.write(f"  broj rolling FD tacaka = {len(fd_roll)}\n")
    f.write(f"  max lag               = {fd_mi_max_lag}\n")
    f.write(f"  bins                  = {fd_mi_bins}\n")
    f.write(f"  max MI                = {fd_max_mi:.8f} bits\n")
    f.write(f"  max MI lag            = {fd_max_mi_lag}\n\n")
    f.write("Shuffled rolling FD max MI referenca:\n")
    f.write(f"  runs                  = {fd_mi_shuffle_runs}\n")
    f.write(f"  mean                  = {shuffle_fd_mi_mean:.8f}\n")
    f.write(f"  std                   = {shuffle_fd_mi_std:.8f}\n")
    f.write(f"  z                     = {shuffle_fd_mi_z:.8f}\n")
    f.write(f"  p(shuffled >= obs)    = {shuffle_fd_mi_p:.8f}\n")
    f.write(f"  interpret.            = {fd_mi_note}\n\n")
    f.write("Top MI lagovi rolling FD:\n")
    f.write(f"  {'lag':<8}{'MI [bits]':>16}\n")
    for lag, val in fd_top_mi_pairs:
        f.write(f"  {lag:<8}{val:>16,.8f}\n")
    f.write("\n")

    elapsed_2c3c = time.time() - T0_2C3C
    f.write(f"Vreme KORAKA 2c3c: {timedelta(seconds=int(elapsed_2c3c))} ({elapsed_2c3c:.1f} s)\n")
    f.write("\nKraj KORAKA 2c3c.\n")

print(f"PNG saved → {PNG_PATH_2C3C}")
print(f"TXT updated → {TXT_PATH}")
print(f"Vreme KORAKA 2c3c: {timedelta(seconds=int(time.time()-T0_2C3C))} "
      f"({time.time()-T0_2C3C:.1f} s)")
print()
"""
KORAK 2c3c (Aparat 2c: Fraktalna dimenzija + Test 3c: Mutual Information)
u 1_KarlWeierstrass_v2.py. Posebna PNG slika, TXT proširen.

KORAK 2c3c — MI nad rolling/local FD nizom,
             MI nad f(t) kao kontrola,
             shuffled max MI referenca.
             Crta rezultat u 1_KarlWeierstrass_v2_2c3c.png.


2c3c: Aparat 2c Fraktalna dimenzija + Test 3c Mutual Information. Radiću MI nad rolling/local FD nizom, MI nad f(t) kao kontrolu i shuffled MI referencu, sa novim PNG _2c3c.

2c3c: MI rolling FD režima, MI f(t) kao kontrolu i shuffled max MI referencu za rolling FD.

Dodato:
1_KarlWeierstrass_v2_2c3c.png
MI nad rolling/local FD nizom
MI nad f(t) kao kontrola
shuffled max MI referenca (500 permutacija)
TXT append i komentar-log blok


KORAK 2c3c: Aparat 2c Fraktalna dimenzija + Test 3c Mutual Information
  rolling FD max MI lag 1..12: 1.512328 bits (lag=8)
  shuffled max MI: mean=1.501248 std=0.106821 z=0.10 p=0.4400
  ⇒ rolling FD nema jak MI signal iznad shuffled reference

PNG saved →   /1_KarlWeierstrass_v2_2c3c.png
TXT updated → /1_KarlWeierstrass_v2.txt
Vreme KORAKA 2c3c: 0:00:15 (15.8 s)
"""



# ─────────────────────────────────────────────────────────────────────
# KORAK 2c3d: Aparat 2c Fraktalna dimenzija + Test 3d Entropy
#   Cilj: izmeriti kompleksnost lokalne fraktalne dimenzije kroz vreme.
#   Rolling FD niz je kratak, zato se sample entropy tretira kao indikativna.
# ─────────────────────────────────────────────────────────────────────
T0_2C3D = time.time()

fd_sampen_m = 2
fd_sampen_r = 0.35
fd_sampen, fd_sampen_a, fd_sampen_b, fd_sampen_n = sample_entropy(
    fd_roll, m=fd_sampen_m, r=fd_sampen_r, max_points=len(fd_roll)
)
fd_control_sampen, _, _, fd_control_sampen_n = sample_entropy(
    lex_idx, m=fd_sampen_m, r=0.2, max_points=1200
)

fd_pe_orders = [3, 4, 5]
fd_pe_rows = []
fd_control_pe_rows = []
for order in fd_pe_orders:
    pe, pe_norm, patterns = permutation_entropy(fd_roll, order=order, delay=1)
    fd_pe_rows.append((order, pe, pe_norm, patterns, math.factorial(order)))
    pe_f, pe_f_norm, patterns_f = permutation_entropy(lex_idx, order=order, delay=1)
    fd_control_pe_rows.append((order, pe_f, pe_f_norm, patterns_f, math.factorial(order)))

rng_2c3d = np.random.default_rng(53)
fd_entropy_shuffle_runs = 500
shuffle_fd_sampen = []
shuffle_fd_pe4 = []
for _ in range(fd_entropy_shuffle_runs):
    shuffled_fd = rng_2c3d.permutation(fd_roll)
    se, _, _, _ = sample_entropy(
        shuffled_fd, m=fd_sampen_m, r=fd_sampen_r, max_points=len(shuffled_fd)
    )
    _, pe4_norm_shuf, _ = permutation_entropy(shuffled_fd, order=4, delay=1)
    if np.isfinite(se):
        shuffle_fd_sampen.append(se)
    shuffle_fd_pe4.append(pe4_norm_shuf)

shuffle_fd_sampen = np.asarray(shuffle_fd_sampen, dtype=float)
shuffle_fd_pe4 = np.asarray(shuffle_fd_pe4, dtype=float)
fd_pe4_norm = fd_pe_rows[1][2]

if len(shuffle_fd_sampen) > 0 and np.isfinite(fd_sampen):
    shuffle_fd_sampen_mean = float(shuffle_fd_sampen.mean())
    shuffle_fd_sampen_std = float(shuffle_fd_sampen.std(ddof=1))
    shuffle_fd_sampen_p_low = float(np.mean(shuffle_fd_sampen <= fd_sampen))
    shuffle_fd_sampen_z = (
        (fd_sampen - shuffle_fd_sampen_mean) / (shuffle_fd_sampen_std + 1e-12)
    )
else:
    shuffle_fd_sampen_mean = float("nan")
    shuffle_fd_sampen_std = float("nan")
    shuffle_fd_sampen_p_low = float("nan")
    shuffle_fd_sampen_z = float("nan")

shuffle_fd_pe4_mean = float(shuffle_fd_pe4.mean())
shuffle_fd_pe4_std = float(shuffle_fd_pe4.std(ddof=1))
shuffle_fd_pe4_p_low = float(np.mean(shuffle_fd_pe4 <= fd_pe4_norm))
shuffle_fd_pe4_z = (fd_pe4_norm - shuffle_fd_pe4_mean) / (shuffle_fd_pe4_std + 1e-12)

if (np.isfinite(shuffle_fd_sampen_p_low) and shuffle_fd_sampen_p_low <= 0.05) or (
    shuffle_fd_pe4_p_low <= 0.05
):
    fd_entropy_note = "entropy rolling FD je niza od shuffled reference (moguca struktura)"
else:
    fd_entropy_note = "entropy rolling FD je blizu shuffled reference"

print()
print("KORAK 2c3d: Aparat 2c Fraktalna dimenzija + Test 3d Sample / Permutation entropy")
print(f"  Sample entropy rolling FD: {fd_sampen:.4f} "
      f"(m={fd_sampen_m}, r={fd_sampen_r}, n={fd_sampen_n})")
print(f"  Sample entropy f(t) kontrola: {fd_control_sampen:.4f} (n={fd_control_sampen_n})")
print(f"  Permutation entropy rolling FD order=4: {fd_pe4_norm:.4f} normalizovano")
print(f"  shuffled SampEn: mean={shuffle_fd_sampen_mean:.4f} "
      f"std={shuffle_fd_sampen_std:.4f} z={shuffle_fd_sampen_z:.2f} "
      f"p_low={shuffle_fd_sampen_p_low:.4f}")
print(f"  shuffled PE4: mean={shuffle_fd_pe4_mean:.4f} std={shuffle_fd_pe4_std:.4f} "
      f"z={shuffle_fd_pe4_z:.2f} p_low={shuffle_fd_pe4_p_low:.4f}")
print(f"  ⇒ {fd_entropy_note}")
print()

fig2c3d, ax2c3d = plt.subplots(1, 3, figsize=(16, 5))
fig2c3d.suptitle("KORAK 2c3d: Fraktalna dimenzija + Sample / Permutation entropy",
                 fontsize=13, fontweight="bold")

orders = np.array([row[0] for row in fd_pe_rows], dtype=int)
fd_pe_norms = np.array([row[2] for row in fd_pe_rows], dtype=float)
fd_control_pe_norms = np.array([row[2] for row in fd_control_pe_rows], dtype=float)
ax2c3d[0].plot(orders, fd_pe_norms, "o-", color="purple", label="rolling FD")
ax2c3d[0].plot(orders, fd_control_pe_norms, "o-", color="steelblue", label="f(t)")
ax2c3d[0].set_ylim(0, 1.05)
ax2c3d[0].set_title("Normalizovana permutation entropy")
ax2c3d[0].set_xlabel("order")
ax2c3d[0].set_ylabel("PE / max PE")
ax2c3d[0].legend(fontsize=8)
ax2c3d[0].grid(True, alpha=0.25)

if len(shuffle_fd_sampen) > 0:
    ax2c3d[1].hist(shuffle_fd_sampen, bins=20, color="lightgray", edgecolor="white")
    ax2c3d[1].axvline(fd_sampen, color="crimson", linewidth=2,
                      label=f"observed={fd_sampen:.3f}")
    ax2c3d[1].axvline(shuffle_fd_sampen_mean, color="black", linestyle="--",
                      label=f"shuffle mean={shuffle_fd_sampen_mean:.3f}")
ax2c3d[1].set_title("Shuffled rolling FD Sample Entropy")
ax2c3d[1].set_xlabel("Sample entropy")
ax2c3d[1].set_ylabel("broj")
ax2c3d[1].legend(fontsize=8)

ax2c3d[2].hist(shuffle_fd_pe4, bins=20, color="lightgray", edgecolor="white")
ax2c3d[2].axvline(fd_pe4_norm, color="crimson", linewidth=2,
                  label=f"observed={fd_pe4_norm:.3f}")
ax2c3d[2].axvline(shuffle_fd_pe4_mean, color="black", linestyle="--",
                  label=f"shuffle mean={shuffle_fd_pe4_mean:.3f}")
ax2c3d[2].set_title("Shuffled rolling FD PE order=4")
ax2c3d[2].set_xlabel("normalizovana PE")
ax2c3d[2].set_ylabel("broj")
ax2c3d[2].legend(fontsize=8)

for a in ax2c3d:
    a.spines["top"].set_visible(False)
    a.spines["right"].set_visible(False)

fig2c3d.tight_layout()
fig2c3d.savefig(PNG_PATH_2C3D, dpi=150, bbox_inches="tight")
plt.show()

with open(TXT_PATH, "a", encoding="utf-8") as f:
    f.write("\n")
    f.write("=" * 60 + "\n")
    f.write("KORAK 2c3d: Aparat 2c Fraktalna dimenzija + Test 3d Sample / Permutation entropy\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"  PNG:                  {PNG_PATH_2C3D}\n\n")
    f.write("Sample entropy rolling/local FD niza:\n")
    f.write(f"  m                     = {fd_sampen_m}\n")
    f.write(f"  r                     = {fd_sampen_r}\n")
    f.write(f"  n used                = {fd_sampen_n}\n")
    f.write(f"  SampEn(rolling FD)    = {fd_sampen:.8f}\n")
    f.write(f"  SampEn(f(t) kontrola) = {fd_control_sampen:.8f}\n")
    f.write(f"  A count               = {fd_sampen_a}\n")
    f.write(f"  B count               = {fd_sampen_b}\n\n")
    f.write("Permutation entropy rolling FD:\n")
    f.write(f"  {'order':<8}{'PE bits':>16}{'PE norm':>16}{'patterns':>14}{'max':>10}\n")
    for order, pe, pe_norm, patterns, max_patterns in fd_pe_rows:
        f.write(f"  {order:<8}{pe:>16,.8f}{pe_norm:>16,.8f}{patterns:>14}{max_patterns:>10}\n")
    f.write("\n")
    f.write("Shuffled entropy referenca:\n")
    f.write(f"  runs                  = {fd_entropy_shuffle_runs}\n")
    f.write(f"  SampEn finite runs    = {len(shuffle_fd_sampen)}\n")
    f.write(f"  SampEn mean           = {shuffle_fd_sampen_mean:.8f}\n")
    f.write(f"  SampEn std            = {shuffle_fd_sampen_std:.8f}\n")
    f.write(f"  SampEn z              = {shuffle_fd_sampen_z:.8f}\n")
    f.write(f"  SampEn p_low          = {shuffle_fd_sampen_p_low:.8f}\n")
    f.write(f"  PE4 mean              = {shuffle_fd_pe4_mean:.8f}\n")
    f.write(f"  PE4 std               = {shuffle_fd_pe4_std:.8f}\n")
    f.write(f"  PE4 z                 = {shuffle_fd_pe4_z:.8f}\n")
    f.write(f"  PE4 p_low             = {shuffle_fd_pe4_p_low:.8f}\n")
    f.write(f"  interpret.            = {fd_entropy_note}\n\n")

    elapsed_2c3d = time.time() - T0_2C3D
    f.write(f"Vreme KORAKA 2c3d: {timedelta(seconds=int(elapsed_2c3d))} ({elapsed_2c3d:.1f} s)\n")
    f.write("\nKraj KORAKA 2c3d.\n")

print(f"PNG saved → {PNG_PATH_2C3D}")
print(f"TXT updated → {TXT_PATH}")
print(f"Vreme KORAKA 2c3d: {timedelta(seconds=int(time.time()-T0_2C3D))} "
      f"({time.time()-T0_2C3D:.1f} s)")
print()
"""
KORAK 2c3d (Aparat 2c: Fraktalna dimenzija + Test 3d: Sample / Permutation entropy)
u 1_KarlWeierstrass_v2.py. Posebna PNG slika, TXT proširen.

KORAK 2c3d — sample entropy nad rolling/local FD nizom,
             permutation entropy nad rolling/local FD nizom,
             shuffled entropy referenca.
             Crta rezultat u 1_KarlWeierstrass_v2_2c3d.png.


2c3d: Aparat 2c Fraktalna dimenzija + Test 3d Sample / Permutation entropy. Radiću entropiju nad rolling/local FD nizom, f(t) kao kontrolu i shuffled entropy referencu, sa novim PNG _2c3d.

2c3d: Sample entropy + permutation entropy nad rolling FD nizom, f(t) kao kontrola i shuffled entropy referenca za rolling FD.

Dodato:
1_KarlWeierstrass_v2_2c3d.png
Sample entropy nad rolling/local FD nizom
Permutation entropy nad rolling/local FD nizom
f(t) kao kontrola
shuffled entropy referenca (500 permutacija)
TXT append i komentar-log blok


KORAK 2c3d: Aparat 2c Fraktalna dimenzija + Test 3d Sample / Permutation entropy
  Sample entropy rolling FD: 1.9459 (m=2, r=0.35, n=31)
  Sample entropy f(t) kontrola: 2.1954 (n=1200)
  Permutation entropy rolling FD order=4: 0.8868 normalizovano
  shuffled SampEn: mean=1.8031 std=0.5088 z=0.28 p_low=0.6797
  shuffled PE4: mean=0.8592 std=0.0382 z=0.72 p_low=0.7900
  ⇒ entropy rolling FD je blizu shuffled reference

PNG saved →   /1_KarlWeierstrass_v2_2c3d.png
TXT updated → /1_KarlWeierstrass_v2.txt
Vreme KORAKA 2c3d: 0:00:17 (17.2 s)
"""



# ─────────────────────────────────────────────────────────────────────
# KORAK 2c3e: Aparat 2c Fraktalna dimenzija + Test 3e NIST baterija
#   Cilj: proveriti binarnu nasumicnost lokalnih FD rezima.
#   Bit = 1 znaci da je lokalni FD iznad medijane rolling FD niza.
# ─────────────────────────────────────────────────────────────────────
T0_2C3E = time.time()

fd_bit_threshold = float(np.median(fd_roll))
fd_bits = (fd_roll > fd_bit_threshold).astype(int)
fd_bits_n = len(fd_bits)
fd_bits_ones = int(fd_bits.sum())
fd_bits_zeros = int(fd_bits_n - fd_bits_ones)

fd_monobit_p, fd_monobit_s = nist_monobit(fd_bits)
fd_runs_p, fd_runs_count, fd_runs_pi = nist_runs(fd_bits)
fd_block_p, fd_block_chi2, fd_block_count = nist_block_frequency(fd_bits, block_size=8)
fd_cusum_p, fd_cusum_z, fd_cusum_walk = nist_cumulative_sums(fd_bits)
fd_apen_p, fd_apen_value, fd_apen_chi2, fd_apen_df = nist_approximate_entropy(fd_bits, m=2)

fd_nist_rows = [
    ("Monobit frequency", fd_monobit_p, fd_monobit_s),
    ("Runs", fd_runs_p, fd_runs_count),
    ("Block frequency", fd_block_p, fd_block_chi2),
    ("Cumulative sums", fd_cusum_p, fd_cusum_z),
    ("Approx entropy", fd_apen_p, fd_apen_value),
]
fd_nist_pass_count = int(sum(p > 0.05 for _, p, _ in fd_nist_rows if np.isfinite(p)))
fd_nist_total = int(sum(np.isfinite(p) for _, p, _ in fd_nist_rows))

fd_f_control_bits = nist_bits_from_series(lex_idx)
fd_f_monobit_p, fd_f_monobit_s = nist_monobit(fd_f_control_bits)
fd_f_runs_p, fd_f_runs_count, fd_f_runs_pi = nist_runs(fd_f_control_bits)
fd_f_block_p, fd_f_block_chi2, fd_f_block_count = nist_block_frequency(
    fd_f_control_bits, block_size=128
)
fd_f_cusum_p, fd_f_cusum_z, _ = nist_cumulative_sums(fd_f_control_bits)
fd_f_apen_p, fd_f_apen_value, fd_f_apen_chi2, fd_f_apen_df = nist_approximate_entropy(
    fd_f_control_bits, m=2
)
fd_f_nist_rows = [
    ("Monobit frequency", fd_f_monobit_p, fd_f_monobit_s),
    ("Runs", fd_f_runs_p, fd_f_runs_count),
    ("Block frequency", fd_f_block_p, fd_f_block_chi2),
    ("Cumulative sums", fd_f_cusum_p, fd_f_cusum_z),
    ("Approx entropy", fd_f_apen_p, fd_f_apen_value),
]
fd_f_nist_pass_count = int(sum(p > 0.05 for _, p, _ in fd_f_nist_rows if np.isfinite(p)))
fd_f_nist_total = int(sum(np.isfinite(p) for _, p, _ in fd_f_nist_rows))

if fd_nist_pass_count == fd_nist_total:
    fd_nist_note = "rolling FD bitovi prolaze sve NIST-style testove"
elif fd_nist_pass_count >= max(1, fd_nist_total - 1):
    fd_nist_note = "rolling FD bitovi uglavnom prolaze, slab signal za proveru"
else:
    fd_nist_note = "rolling FD bitovi padaju vise NIST-style testova"

print()
print("KORAK 2c3e: Aparat 2c Fraktalna dimenzija + Test 3e NIST baterija")
print(f"  rolling FD bits: n={fd_bits_n}  zeros={fd_bits_zeros}  ones={fd_bits_ones} "
      f"threshold={fd_bit_threshold:.6f}")
for name, p, stat_val in fd_nist_rows:
    print(f"  {name:<20} p={p:.6f}  stat={stat_val}")
print(f"  prolaz rolling FD: {fd_nist_pass_count}/{fd_nist_total}  ⇒ {fd_nist_note}")
print(f"  kontrola f(t) prolaz: {fd_f_nist_pass_count}/{fd_f_nist_total}")
print()

fig2c3e, ax2c3e = plt.subplots(1, 3, figsize=(16, 5))
fig2c3e.suptitle("KORAK 2c3e: Fraktalna dimenzija + NIST-style testovi",
                 fontsize=13, fontweight="bold")

ax2c3e[0].bar(["FD<=median", "FD>median"], [fd_bits_zeros, fd_bits_ones],
              color=["steelblue", "purple"])
ax2c3e[0].set_title("Binarizovani rolling FD rezimi")
ax2c3e[0].set_xlabel("bit")
ax2c3e[0].set_ylabel("broj")
ax2c3e[0].grid(True, alpha=0.2, axis="y")

names = [row[0] for row in fd_nist_rows]
fd_pvals = np.array([row[1] for row in fd_nist_rows], dtype=float)
colors = ["seagreen" if p > 0.05 else "crimson" for p in fd_pvals]
ax2c3e[1].barh(names, fd_pvals, color=colors)
ax2c3e[1].axvline(0.05, color="black", linestyle="--", linewidth=1.2)
ax2c3e[1].set_xlim(0, 1)
ax2c3e[1].set_title("Rolling FD NIST-style p-vrednosti")
ax2c3e[1].set_xlabel("p-value")
ax2c3e[1].grid(True, alpha=0.2, axis="x")

ax2c3e[2].plot(np.arange(1, fd_bits_n + 1), fd_cusum_walk,
               linewidth=1.2, marker="o", markersize=3, color="purple")
ax2c3e[2].axhline(0, color="black", linewidth=0.6)
ax2c3e[2].set_title(f"Rolling FD cumulative sums (z={fd_cusum_z})")
ax2c3e[2].set_xlabel("rolling window index")
ax2c3e[2].set_ylabel("cum sum")
ax2c3e[2].grid(True, alpha=0.25)

for a in ax2c3e:
    a.spines["top"].set_visible(False)
    a.spines["right"].set_visible(False)

fig2c3e.tight_layout()
fig2c3e.savefig(PNG_PATH_2C3E, dpi=150, bbox_inches="tight")
plt.show()

with open(TXT_PATH, "a", encoding="utf-8") as f:
    f.write("\n")
    f.write("=" * 60 + "\n")
    f.write("KORAK 2c3e: Aparat 2c Fraktalna dimenzija + Test 3e NIST baterija\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"  PNG:                  {PNG_PATH_2C3E}\n\n")
    f.write("Binarizacija rolling/local FD niza:\n")
    f.write("  bit = 1 ako je lokalni FD > medijana rolling FD niza, inace 0\n")
    f.write(f"  threshold             = {fd_bit_threshold:.8f}\n")
    f.write(f"  n bits                = {fd_bits_n}\n")
    f.write(f"  zeros                 = {fd_bits_zeros}\n")
    f.write(f"  ones                  = {fd_bits_ones}\n\n")
    f.write("NIST-style testovi nad rolling FD bitovima (prolaz ako p > 0.05):\n")
    f.write(f"  {'test':<22}{'p-value':>14}{'stat':>18}{'pass':>10}\n")
    for name, p, stat_val in fd_nist_rows:
        f.write(f"  {name:<22}{p:>14,.8f}{float(stat_val):>18,.8f}{str(p > 0.05):>10}\n")
    f.write("\n")
    f.write("Detalji rolling FD:\n")
    f.write(f"  Runs pi               = {fd_runs_pi:.8f}\n")
    f.write(f"  Block size            = 8\n")
    f.write(f"  Block count           = {fd_block_count}\n")
    f.write(f"  Approx entropy m      = 2\n")
    f.write(f"  Approx entropy chi2   = {fd_apen_chi2:.8f}\n")
    f.write(f"  Approx entropy df     = {fd_apen_df}\n")
    f.write(f"  pass count            = {fd_nist_pass_count}/{fd_nist_total}\n")
    f.write(f"  interpret.            = {fd_nist_note}\n\n")
    f.write("Kontrola: NIST-style testovi nad f(t) binarizacijom:\n")
    f.write(f"  {'test':<22}{'p-value':>14}{'stat':>18}{'pass':>10}\n")
    for name, p, stat_val in fd_f_nist_rows:
        f.write(f"  {name:<22}{p:>14,.8f}{float(stat_val):>18,.8f}{str(p > 0.05):>10}\n")
    f.write(f"  pass count            = {fd_f_nist_pass_count}/{fd_f_nist_total}\n")
    f.write(f"  f(t) runs pi          = {fd_f_runs_pi:.8f}\n")
    f.write(f"  f(t) block count      = {fd_f_block_count}\n")
    f.write(f"  f(t) approx chi2      = {fd_f_apen_chi2:.8f}\n")
    f.write(f"  f(t) approx df        = {fd_f_apen_df}\n\n")

    elapsed_2c3e = time.time() - T0_2C3E
    f.write(f"Vreme KORAKA 2c3e: {timedelta(seconds=int(elapsed_2c3e))} ({elapsed_2c3e:.1f} s)\n")
    f.write("\nKraj KORAKA 2c3e.\n")

print(f"PNG saved → {PNG_PATH_2C3E}")
print(f"TXT updated → {TXT_PATH}")
print(f"Vreme KORAKA 2c3e: {timedelta(seconds=int(time.time()-T0_2C3E))} "
      f"({time.time()-T0_2C3E:.1f} s)")
print()
"""
KORAK 2c3e (Aparat 2c: Fraktalna dimenzija + Test 3e: NIST baterija)
u 1_KarlWeierstrass_v2.py. Posebna PNG slika, TXT proširen.

KORAK 2c3e — NIST-style testovi nad rolling FD > median(FD) bitovima:
             monobit, runs, block frequency, cumulative sums,
             approximate entropy; f(t) binarizacija kao kontrola.
             Crta rezultat u 1_KarlWeierstrass_v2_2c3e.png.


2c3e: Aparat 2c Fraktalna dimenzija + Test 3e NIST baterija. Radiću NIST-style testove nad binarizovanim rolling/local FD nizom, uz f(t) kontrolu, novi PNG _2c3e.

2c3e: NIST-style testovi nad rolling FD > median(FD) bitovima, plus kontrolni NIST nad f(t) binarizacijom.

Dodato:
1_KarlWeierstrass_v2_2c3e.png
NIST-style testovi nad rolling FD > median(FD) bitovima
f(t) NIST binarizacija kao kontrola
Monobit, Runs, Block frequency, Cumulative sums, Approximate entropy
TXT append i komentar-log blok


KORAK 2c3e: Aparat 2c Fraktalna dimenzija + Test 3e NIST baterija
  rolling FD bits: n=31  zeros=16  ones=15 threshold=1.997626
  Monobit frequency    p=0.857462  stat=0.1796053020267749
  Runs                 p=0.593634  stat=14
  Block frequency      p=0.682270  stat=1.5
  Cumulative sums      p=0.590014  stat=3
  Approx entropy       p=0.636552  stat=0.6785765664660117
  prolaz rolling FD: 5/5  ⇒ rolling FD bitovi prolaze sve NIST-style testove
  kontrola f(t) prolaz: 4/5

PNG saved →   /1_KarlWeierstrass_v2_2c3e.png
TXT updated → /1_KarlWeierstrass_v2.txt
Vreme KORAKA 2c3e: 0:00:15 (15.5 s)
"""



"""
============================================================
ANALIZA
============================================================

Napomena:
  POGODNO znaci: rezultat daje merljiv signal koji moze da udje u sledeci
  prediktivni model kao feature / filter / tezina.
  NIJE POGODNO znaci: test je vise dijagnosticki ili je signal preslab,
  nestabilan, trivijalan ili previse blizu shuffled/random reference.

Tabela 15 modela (aparati 2a/2b/2c x testovi 3a/3b/3c/3d/3e):

| Model | Aparat | Test | Prednosti | Mane | Analiza | Zakljucak |
|---|---|---|---|---|---|---|
| 2a3a | Brownovo kretanje | Hurst | Brown-putanja ima H=0.6294 i dobar R2=0.9950 | Shuffled Brown H je mnogo visi; H(dX)=0.0690 slab fit | Nije potvrda klasicnog Browna; vise pokazuje da putanja nije standardna Brown referenca | NIJE POGODNO |
| 2a3b | Brownovo kretanje | ACF | Jak lag-1 signal: max abs ACF=0.5019; Ljung-Box p=0; shuffled p=0 | Signal je dominantno lokalni i linearan; moze biti posledica lex-rank geometrije | Najjaci Brown test za praktican feature: prethodni inkrement/lag-1 nosi informaciju | POGODNO |
| 2a3c | Brownovo kretanje | Mutual information | Jak MI na lag=1: 0.2637 bits; shuffled p=0 | MI opada posle lag-1; zavisi od binovanja | Potvrdjuje nelinearnu/linearno-neuhvatljivu vezu inkremenata; koristan signal za model | POGODNO |
| 2a3d | Brownovo kretanje | Sample / permutation entropy | PE4 je znatno niza od shuffled reference; p_low=0 | SampEn nije znacajno nizi; deo signala dolazi iz ordinalnog obrasca | Entropija pokazuje red u redosledu inkremenata; pogodno kao dodatni skor/feature | POGODNO |
| 2a3e | Brownovo kretanje | NIST baterija | Runs i ApproxEntropy padaju; 3/5 prolaz | Monobit i block prolaze; test radi nad binarnim znakom, gubi amplitudu | Binarni niz inkremenata nije potpuno random; koristan za detekciju rezima, ne samostalno | POGODNO |
| 2b3a | Hurst/R-S | Hurst | H(f)=0.5931; rolling H mean=0.5974; 30/31 prozora H>0.55; shuffled p_high=0.03 | Efekat je umeren; R/S moze preceniti memoriju kod konacnih nizova | Najbolji dokaz perzistentnog rezima u f(t); dobar kandidat za trend/regime feature | POGODNO |
| 2b3b | Hurst/R-S | ACF | Rolling H ima ACF signal; Ljung-Box p=0.0154; shuffled p=0.028 | Samo 31 rolling tacaka; mali uzorak za ACF | Lokalni H rezimi se grupisu kroz vreme; korisno za rezimsku predikciju | POGODNO |
| 2b3c | Hurst/R-S | Mutual information | MI rolling H je visok nominalno | Shuffled referenca je jos visa; p=0.922 | MI rolling H nije iznad randomizovane reference; nema dodatnu prediktivnu vrednost | NIJE POGODNO |
| 2b3d | Hurst/R-S | Sample / permutation entropy | Entropije su izmerene i stabilne | Shuffled reference ne potvrdjuju odmak; PE4 p_low=0.588 | Kompleksnost rolling H nije posebna u odnosu na shuffled; dijagnosticki, ne prediktivno | NIJE POGODNO |
| 2b3e | Hurst/R-S | NIST baterija | Rolling H bitovi padaju 5/5 testova; svi prozori H>0.5 | Signal je skoro trivijalan jer zeros=0, ones=31 | Jako pokazuje perzistentan rezim, ali kao binarni test je pregrub; ipak koristan kao rezimski filter | POGODNO |
| 2c3a | Fraktalna dimenzija | Hurst | Rolling FD uporedjen sa 2-H; corr(FD,H)=-0.4357 | Globalni FD=1.9988 je skoro shuffled; p_high=0.85 | FD je ekstremno hrapav skoro kao random; H veza ne daje stabilan prediktivni signal | NIJE POGODNO |
| 2c3b | Fraktalna dimenzija | ACF | Rolling FD ima ACF signal; Ljung-Box p=0.000724; shuffled p=0.014 | FD vrednosti su veoma zbijene oko ~2; mala dinamika | Lokalna hrapavost ima memoriju; moze biti pomocni rezimski feature | POGODNO |
| 2c3c | Fraktalna dimenzija | Mutual information | MI rolling FD postoji nominalno | Shuffled p=0.44; nije iznad reference | Nema dokaz dodatne nelinearne informacije u FD rezimu | NIJE POGODNO |
| 2c3d | Fraktalna dimenzija | Sample / permutation entropy | Entropije su uredno izmerene | Shuffled p vrednosti nisu znacajne; PE4 p_low=0.79 | Entropija rolling FD je blizu shuffled reference; nema prediktivni signal | NIJE POGODNO |
| 2c3e | Fraktalna dimenzija | NIST baterija | Rolling FD bitovi prolaze 5/5 testova | Prolazak NIST-a znaci nasumicnost, ne strukturu | FD binarni rezim izgleda random; nije dobar za predikciju sledece kombinacije | NIJE POGODNO |

Kratak zakljucak:
  POGODNO za sledecu fazu predikcije (8/15):
    1) 2a3b - Brown inkrementi + ACF
    2) 2a3c - Brown inkrementi + Mutual information
    3) 2a3d - Brown inkrementi + sample/permutation entropy
    4) 2a3e - Brown inkrementi + NIST baterija
    5) 2b3a - Hurst/R-S globalni i rolling H
    6) 2b3b - ACF rolling H rezima
    7) 2b3e - NIST rolling H rezima
    8) 2c3b - rolling FD ACF kao pomocni rezimski signal

  NIJE POGODNO za predikciju sledece kombinacije (7/15):
    1) 2a3a - Brown + Hurst: ne potvrdjuje standardni Brown signal
    2) 2b3c - Hurst + MI: nije iznad shuffled reference
    3) 2b3d - Hurst + entropy: blizu shuffled reference
    4) 2c3a - FD + Hurst: globalni FD skoro kao shuffled
    5) 2c3c - FD + MI: nije iznad shuffled reference
    6) 2c3d - FD + entropy: blizu shuffled reference
    7) 2c3e - FD + NIST: prolazi NIST, tj. izgleda random

  Operativno:
    Za sledeci prediktivni model ne treba koristiti sve jednako.
    Osnovu treba graditi iz Brown inkremenata (lag/ACF/MI) + Hurst rezima,
    a fraktalnu dimenziju koristiti samo tamo gde pokazuje lokalnu ACF memoriju.
"""

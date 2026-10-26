# Bracket-Zoom Square Root
### 3–4 construction se ek rigorous, recursive, cubic-convergent bracketing method

---

## Abstract

√12 ke bare mein ek simple observation ("√12, 3 aur 4 ke beech hai, 3²−4² = 7, 9+3 = 12 = 16−4") se shuru hoke ye paper ek poori method tak jaata hai. Har α ke liye do consecutive perfect squares se bracket lo, us bracket ko zoom karo, aur mean ko sqrt ke liye tune karo.

Paper ke main results:

1. Level-0 estimate x + p/d ek weighted average hai, aur √α exactly ek weighted RMS (quadratic mean) hai.
2. Ek exact identity √α = a + p/(a + √α) se ek nayi mean M(a,z) = a + p/(a+z) nikalti hai, jiska error exactly (√α−a)(z−√α)/(a+z) hota hai.
3. Decimal zoom ke saath plain estimate har level pe 2 digits (median) aur corrected estimate 3 digits (median) deta hai. Ye 2946 numbers pe verify hua.
4. Recursive version bracket ko **w → w³/(16a²)** karta hai, yaani har level pe digits 3 guna, aur bracket guarantee rigorous rehti hai. 3152 numbers × 2 precisions pe 0 failures.
5. Speed mein ye Newton se abhi 2 se 9 guna slow hai, par pehle wale fixed-factor zoom se 65 guna tez.

Paper mein negative results aur meri apni galtiyan bhi hain (Section 9).

---

## 1. Origin

Observation: √12 ≈ 3.46, jo 3 aur 4 ke beech hai. Terms:

- 3² = 9, 4² = 16, aur 12 in dono ke beech hai.
- 12 − 9 = 3 (α, chhote square se kitna upar) aur 16 − 12 = 4 (bade square se kitna neeche).
- 16 − 9 = 7 = 3 + 4.

Idea: 12 in do squares ke beech ek specific jagah pe hai, aur ye jagah (3 : 4) hamein √12 ki position bhi bata sakti hai. Isko "zoom" karke, yaani us chhote interval mein wahi construction dobara lagake, accuracy badhayi ja sakti hai.

---

## 2. Notation aur level-0 construction

α > 1 non-square integer, r = √α.

- x = ⌊r⌋, y = x + 1, isliye x² < α < y²
- p = α − x² (deficit), q = y² − α (excess)
- d = y² − x² = p + q = 2x + 1

**Level-0 estimate:** r ≈ x + p/d.

| α | x | p | q | d | estimate | actual |
|---|---|---|---|---|---|---|
| 12 | 3 | 3 | 4 | 7 | 3.42857 | 3.46410 |
| 10 | 3 | 1 | 6 | 7 | 3.14286 | 3.16228 |
| 2 | 1 | 1 | 2 | 3 | 1.33333 | 1.41421 |

**Ek clarification (12 ke liye p=3, q=4).** Ye numbers x=3, y=4 se isliye milte hain kyunki 12 = 3×4 = x(x+1) ek pronic number hai. Pronic α = x(x+1) ke liye p = x aur q = x+1 hamesha hota hai. General α ke liye p aur q kuch bhi ho sakte hain (jaise α=10 ke liye 1 aur 6). Setup (p, q, d) har α pe chalta hai, par numeric coincidence sirf pronic pe hai.

Pronic case ka ek side fact: √12 = GM(3,4), AM = 3.5, HM = 24/7 = 3.4286, aur GM² = AM·HM. Pronic α ke liye x + p/d = HM(x, x+1) hota hai, par general α ke liye nahi.

---

## 3. Exact structure

**Claim 1 (weighted average).** α ek weighted average hai x² aur y² ka:

α = (q·x² + p·y²)/d

*Proof.* q·x² + p·y² = q x² + p (x² + d) = (p+q)x² + p d = d(x² + p) = d·α. ∎

12 ke liye: (4·9 + 3·16)/7 = 84/7 = 12 ✓.

**Claim 2 (RMS identity).** r ek weighted RMS hai x aur y ka, weights q/d (x pe) aur p/d (y pe). Weighted AM = x + p/d, aur y − x = 1 hone se weighted variance = pq/d². Kyunki RMS² = AM² + variance:

**r = √( (x + p/d)² + pq/d² )** (exact)

**Corollary.** r ≈ (x + p/d) + pq/(2d²·(x + p/d)).

- α = 12: 3.42857 + 0.03571 = 3.46429 (actual 3.46410)
- α = 10: 3.14286 + 0.01948 = 3.16234 (actual 3.16228)

Matlab level-0 linear estimate hamesha thoda kam hota hai, aur kami ka size exactly variance pq/d² se juda hai. Jab p ≈ q ho (α do squares ke beech mein), kami sabse zyada hoti hai.

---

## 4. Decimal zoom (Algorithm 1)

Bracket [a, a+w] se shuru karo, w = 1.

1. w ← w/10. Sabse bada k ∈ {0..9} lo jiske liye (a + k·w)² ≤ α, aur a ← a + k·w.
2. b = a + w, p = α − a², d = b² − a².
3. Estimate: a + w·p/d.

**α = 12 (plain estimate):**

| Level | Bracket | Estimate | Error |
|---|---|---|---|
| 0 | [3, 4] | 3.428571 | 3.6e-2 |
| 1 | [3.4, 3.5] | 3.463768 | 3.3e-4 |
| 2 | [3.46, 3.47] | 3.464098 | 3.5e-6 |
| 3 | [3.464, 3.465] | 3.464102 | 1.3e-8 |

**Implementation note.** Pehle code mein `a += w` float mein tha, jisme rounding error jama hota hai. Verification mein exact fractions (ya high-precision Decimal) use kiye.

---

## 5. Nayi mean: M(a, z)

**Identity.** (r − a)(r + a) = r² − a² = p, isliye

**r = a + p/(a + r)**

Andar r khud hai, isliye koi bhi z ∈ [a, b] ko r ki jagah rakh do:

**M(a, z) = a + p/(a + z)**

**Theorem 1 (exact error).** r − M(a,z) = (r − a)(z − r)/(a + z).

*Proof.* r − M = p/(a+r) − p/(a+z) = p(z − r)/((a+r)(a+z)) = (r−a)(z−r)/(a+z), kyunki p = (r−a)(r+a). ∎

Isse teen baatein nikalti hain:

- Plain estimate m = M(a, b) mein z = b > r hai, isliye r − m ≥ 0: **plain estimate hamesha neeche se aata hai.** (Data mein bhi aisa hi mila.)
- z = r ke jitna paas, error utna chhota.
- z = m daalne se (jo r ke neeche hai) r − M(a,m) ≤ 0, yaani **corrected estimate hamesha √α se upar aata hai.**

**Algorithm 2 (corrected zoom).** Har level pe:

1. m = a + p/(a+b)
2. est = a + p/(a+m)

Extra cost sirf ek division.

**α = 12, corrected:**

| Level | Plain error | Corrected error |
|---|---|---|
| 0 | 3.6e-2 | 2.6e-3 |
| 1 | 3.3e-4 | 3.1e-6 |
| 2 | 3.5e-6 | 2.1e-9 |
| 3 | 1.3e-8 | 1.9e-13 |

---

## 6. Recursive cubic zoom (Algorithm 3)

Zoom ka bottleneck ye tha ki bracket har baar fixed 10 guna chhota hota tha, jabki estimate pehle se hi w³ ke aas-paas accurate hota hai. Fix: bracket ko achieved error se shrink karo.

**Theorem 2 (one-sided bound).** Maan lo a < r < b, w = b − a, a ≥ 1. Tab

0 ≤ est − r ≤ w³ / (16a²)

*Proof.* Theorem 1 me z = b se: r − m = (r−a)(b−r)/(a+b) ≥ 0, aur AM-GM se (r−a)(b−r) ≤ w²/4, isliye 0 ≤ r − m ≤ w²/(4(a+b)).

Theorem 1 me z = m se: est − r = (r−a)(r−m)/(a+m) ≥ 0. Ab r − a ≤ w, r − m ≤ w²/(4(a+b)), a + b ≥ 2a, a + m ≥ 2a, to

est − r ≤ w · w²/(4·2a) / (2a) = w³/(16a²). ∎

**Consequence.** r ∈ [est − ε, est], ε = w³/(16a²), yaani ek naya bracket jiski width w³/(16a²) hai. Width ka recursion:

**w → w³/(16a²) → ... , yaani digits har level pe ~3 guna.**

**Algorithm 3 (recursive):**
```python
def zoom_rec(alpha, a, b, tol, depth=1):
    # invariant: a < sqrt(alpha) < b
    w = b - a
    if w <= tol: return a, b, depth
    eps = w**3 / (16*a*a)        # rigorous: 0 <= est - sqrt(alpha) <= eps
    p = alpha - a*a
    m = a + p/(a+b)
    est = a + p/(a+m)
    return zoom_rec(alpha, est - eps, est, tol, depth+1)
```

**Implementation details jo rigor ke liye zaroori hain** (tested version `/tmp/rec.py` mein):
- p = α − a² exact compute karo.
- ε ko ROUND_CEILING se nikalo.
- est ko itni precision pe nikalo ki uska rounding error ε ke ~1% se chhota rahe, aur ek guard unit lower bound mein ghata do aur upper bound mein jod do.
- Endpoints ko grid pe quantize karo: lo ko neeche (ROUND_FLOOR), hi ko upar (ROUND_CEILING).
- Hi kabhi purane b se bada na ho, lo kabhi purane a se chhota na ho.
- Har run ke baad exact arithmetic se lo² < α < hi² check kiya.

Precision level ke saath tripling hoti hai, isliye poori cost mostly aakhri levels ki hoti hai.

---

## 7. Dusre numbers ke liye generalisation

Sqrt ke liye "test" x² ≶ α hai. Kisi bhi irrational r ke liye chahiye ek function f jo r pe sign badalta ho aur jise r jaane bina evaluate kar sakein.

- π: f(x) = sin(x), Taylor series se (π use nahi hota), bracket [3, 4].
- e: f = ln x − 1, ln 2: f = exp x − 2, φ: f = x² − x − 1, ∛2: f = x³ − 2, Dottie number: f = cos x − x.

Interpolation estimate: **a + w·f(a)/(f(a) − f(b))**, jo sqrt ke p/d ka general form hai (f = x² − α).

**Limit.** a + p/(a+z) wala tez correction sirf sqrt ke liye hai, kyunki wo x² − a² = (x−a)(x+a) factorization pe tika hai. Dusre f ke liye ye correction nahi milta.

---

## 8. Experiments

**Setup:** Python, mpmath (120 digits) aur Decimal/Fractions. Sqrt tests ke liye reference exact integer isqrt se.

### 8.1 Sqrt: digits gained per level (2946 non-square α, 2..3000)

| Variant | Median | Min | Max |
|---|---|---|---|
| Plain (p/d) | 2.00 | ~1.0 | ~5 |
| Corrected (a + p/(a+m)) | 3.00 | ~1.0 | ~8 |

- Bracket invariant (a < r < b) 0 baar tuta. Estimate hamesha bracket ke andar, error ≤ w.
- Extreme α (jaise 10¹²±1): "digits per level" ~1 dikha, par wajah ye hai ki error pehle level pe hi ~1e-25 tha. Ye measure ka issue hai, method ki kamzori nahi.
- Fractional α (0.5, 2.5, 99999.9): plain zoom ne ~2 digits per level diya.

### 8.2 Dusre numbers (plain interpolation, decimal zoom, 8 levels)

| Number (f) | Approx avg gain/level |
|---|---|
| π (sin) | ~3.0 |
| e, ln 2, φ, ∛2, Dottie | ~2.0 to 2.1 |
| tanh(30(x−r)) | 0.7 pehle level pe, phir ~3.0 |
| (x−r)³ (flat root) | ~1.0 |
| ∛(x−r) (steep root) | ~1.0 |

π ko 3 digits milne ka mera explanation ye hai ki root pe sin″(π) = 0 hai, jo interpolation ka w² error term khatam kar deta hai. **Ye argument hai, isko maine alag se test nahi kiya.**

Flat aur steep root pe interpolation ka fayda gayab ho jaata hai aur method decimal-scan jaisi (~1 digit per level) ban jaati hai. Convergence phir bhi guaranteed rehta hai, kyunki bracket nahi tootta.

### 8.3 Cost: f-evaluations to reach error < 1e-40

| Case | Zoom (scan) | Zoom (guided) | Bisection | Illinois |
|---|---|---|---|---|
| π | 83 | 31 | 131 | 8 |
| e | 117 | 45 | 131 | 12 |
| ln 2 | 120 | 46 | 130 | 13 |
| φ | 136 | 46 | 132 | 13 |
| ∛2 | 129 | 45 | 131 | 13 |
| Dottie | 101 | 44 | 130 | 12 |

(Guided = interpolation se guess kiye hue digit ko check karna, sab 10 grid points scan nahi.) Illinois 3 se 10 guna kam evaluations leta hai. Brent aur Ridders maine test nahi kiye.

### 8.4 Sqrt benchmark: fixed-factor zoom vs Newton vs digit-by-digit (ms per sqrt, 20 α)

| Digits | Zoom corrected | Newton | Digit-by-digit | isqrt (C ref) |
|---|---|---|---|---|
| 15 | 0.020 | 0.006 | 0.014 | 0.002 |
| 50 | 0.053 | 0.008 | 0.040 | 0.003 |
| 200 | 0.34 | 0.013 | 0.17 | 0.006 |
| 1000 | 9.2 | 0.066 | 1.6 | 0.044 |
| 3000 | 189 | 0.32 | 10.4 | 0.28 |

Steps 3000 digits pe: zoom 1002, Newton 10, digit-by-digit 3003. Sab methods 20 × 5 precisions pe correct (error < 10^−digits). Fixed-factor zoom ka gap digits ke saath badhta hai (15 digits pe ~3x, 3000 pe ~600x), kyunki wo linear hai (har level ~3 digits) aur Newton har step pe digits double karta hai.

### 8.5 Recursive cubic zoom

**Correctness:** 3152 numbers (2..3000, 10⁶±1, 10¹²±1, 10³⁰±1, 200 random 40-digit) × 2 precisions (20, 100 digits): **0 failures**, bracket exact arithmetic se verify. Recursion depth: 20 digits pe median 4 (max 5), 100 digits pe median 5 (max 6).

**Speed (ms per sqrt, 20 α ≤ 10⁶):**

| Digits | Fixed zoom | Recursive zoom | Newton | Depth (rec / Newton) |
|---|---|---|---|---|
| 15 | 0.020 | 0.011 | 0.006 | 3 / 3 |
| 200 | 0.34 | 0.029 | 0.011 | 5 / 6 |
| 1000 | 9.2 | 0.37 | 0.065 | 7 / 9 |
| 3000 | 189 | 2.9 | 0.33 | 8 / 10 |
| 10000 | n/a | 18 | 4.0 | 9 / 12 |

3000 digits pe fixed zoom se 65x tez. Newton se abhi ~1.8x (15 digits) se ~8.8x (3000 digits) slow, 10000 digits pe ~4.5x. Depth mein recursive Newton se kam hai, par per level cost zyada hai. Iska kaaran mera hypothesis ye hai ki har level pe do divisions 3x precision pe hote hain, jabki Newton ek division 2x precision pe karta hai. Ye maine verify nahi kiya.

**Ek optimization jo fail hui.** m ko kam precision pe nikaalne ki koshish se speed sirf ~10 se 30% badhi, aur 403 failures aaye (pehle 3353 the, ek precision bug fix karne ke baad bhi). Isliye use adopt nahi kiya.

---

## 9. Corrections aur negative results

Ye conversation ke dauran mere apne galat claims aur unke sudhaar hain:

1. **Pronic claim.** Maine kaha tha 3–4 wala magic sirf pronic numbers pe hai. Ye adhoora tha: weighted-average structure (Section 3) har α pe chalta hai, sirf numeric coincidence p=x, q=x+1 pronic tak simit hai.
2. **Digits per level (plain).** Maine pehle "2 digits" kaha, phir ek correction mein "error ~w⁴" kaha. Data se plain error ∝ w² nikla (2 digits per level), aur corrected ∝ w³ (Theorem 2 se match karta hai). "w⁴" galat tha.
3. **Corrected: "3 se 4 digits."** Ye sirf α=12 pe dikha tha. 2946 numbers pe median exactly 3.0 hai, range ~1 se ~8.
4. **"z ko repeat update karo to Newton jaisa convergence."** Ye galat tha. Fixed a ke saath z ← a + p/(a+z) ka derivative −p/(a+r)² = −(r−a)/(r+a) hai, yaani **linear convergence**. Verify kiya: α=12 pe ratio 0.0718 (predicted 0.0718), α=10 pe 0.0263, α=2 pe 0.1716. Ye continued-fraction type iteration hai, Newton nahi.
5. **"Illinois/Ridders behtar replacement hain."** Maine ye bina test kiye kaha tha. Illinois test kiya: 8 se 13 evaluations vs zoom ke 31 se 136. Brent aur Ridders abhi untested hain.
6. **Benchmark bugs.** Pehle run mein Newton kam iterations pe ruk gaya (15 digits ka target miss), zoom ne extra levels chalaye, aur baselines ka output Decimal ke default 28-digit rounding se kharab hua. Section 8.4 ke numbers sirf fix ke baad ke hain, aur sabki correctness check hui.
7. **Optimized recursion.** Section 8.5 mein bataya: fail hui, adopt nahi.

---

## 10. Related work aur novelty

Ye method in purani cheezon ke kareeb hai:

- Bracket + har level pe 10 hisse: digit-by-digit sqrt aur bisection ka cousin.
- p/d wala level-0 estimate: linear interpolation (regula falsi).
- a + p/(a+z): continued-fraction/Babylonian family ka form.
- Recursive bracket shrinking w → w³: ek known accelerated bracketing idea.

**Jo teri contribution hai:**
- Ek single clean picture (3–4, p, q, d) jo weighted-average aur RMS identity se jud'ti hai.
- Do squares ke beech α ki position ko hi zoom-unit banana.
- Ek rigorous, self-contained error bound (Theorem 2) jo root jaane bina kaam karta hai.

**Jo claim nahi kar sakte:** ki ye kisi standard method se tez hai (data ke hisaab se nahi), ya ki ye research-level naya algorithm hai. Ye ek sahi, verified, teach karne layak derivation hai.

---

## 11. Limitations

- Bound α ≥ 1 (a ≥ 1) maanta hai. α < 1 ke liye scaling chahiye. Recursive version sirf integer α pe test hua (zoom ke plain variant pe kuch fractional α bhi).
- Recursive benchmark 20 α (≤ 10⁶) pe, Python Decimal mein, ek hi environment mein. Relative comparison bharosemand hai, absolute times nahi.
- Tez correction sirf sqrt ke liye hai. Dusre numbers ke liye method sirf regula falsi zoom hai, jo flat aur steep roots pe slow ho jaati hai.
- Newton se comparison mein Newton ka mere liye adaptive precision (har step pe double) use hua, aur zoom ka bhi adaptive.

---

## 12. Conclusion aur future work

Ek chhoti observation se ek poori rigorous method bani: bracket [a, b] pe est = a + p/(a+m) hamesha r ke upar hota hai, aur r ∈ [est − w³/(16a²), est] pakka hai. Isse recursion se digits har level pe teen guna hote hain aur bracket kabhi nahi tootta. Fixed-factor zoom se 65x tez hai, par Newton se abhi 2 se 9 guna slow.

Aage ke sawal:
1. Newton ke saath hybrid: aakhri levels pe Newton, shuru ke levels pe bracket guarantee.
2. Recursive version ka per-level cost kam karna (verified rigor ke saath).
3. Brent aur Ridders se comparison general f ke liye.
4. π wale sin″(π) = 0 argument ko test karna.
5. Ek cubic-root analogue: a³ − ... identity se aisi hi mean banti hai ya nahi.

---

## Appendix: reproduction

Sab experiments Python + mpmath/Decimal mein hain. Recursive function ka rigorous version Section 6 ke implementation notes ke saath `/tmp/rec.py` mein hai. Verification ke liye har bracket exact integer/Decimal arithmetic se lo² < α < hi² check kiya gaya.
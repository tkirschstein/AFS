---
title: "MILP Formulation for Agroforestry Biomass Supply Chain (v10)"
author: "SmartAgroforst"
date: "January 2026"
output:
  html_document:
    df_print: paged
  pdf_document: default
documentclass: article
classoption:
- 11pt
- a4paper
geometry: margin=1.25in
fontfamily: libertine
fontsize: 11pt
linestretch: 1.2
numbersections: true
toc: true
toc-depth: 3
colorlinks: true
urlcolor: blue
linkcolor: darkblue
subtitle: "Multi-Cycle Harvesting with Age Lag Constraints"
---


# 1. Index Sets

\begin{align}
i \in \mathcal{I} &= \{1, \ldots, n_s\} && \text{Agroforestry sites} \\
j \in \mathcal{J} &= \{1, \ldots, n_j\} && \text{Storage/processing facilities} \\
k \in \mathcal{K} &= \{1, \ldots, n_k\} && \text{Consumer sites} \\
t \in \mathcal{T} &= \{1, \ldots, T_{\max}\} && \text{Time periods (years)} \\
t \in \mathcal{T^+} &= \{0, \ldots, T_{\max}+1\} && \text{Time periods (years)} \\
t \in \mathcal{T^{harv}} &= \{A^{min}+1, \ldots, T_{\max}\} && \text{ Harvesting periods (years)} \\
p \in \mathcal{P} &= \{1, 2, 3\} && \text{Product types (chemical, pulp, energy)} \\
(s,t) \in \mathcal{S} &= \left\{(s,t) | s,t \in 0,...,(T_{\max}+1) \wedge (t-s) \in \left\{A^{min},...,A^{max}\right\} \right\} && \text{Consecutive harvesting}
\end{align}

---

# 2. Decision Variables

## 2.1 Binary Variables

\begin{align}
z_{ist} & \in \{0, 1\} && \text{Indicator if at site } i \text{ in periods} t,s \text{ are consecutive harvests with } (s,t) \in \mathcal{S}
\end{align}

## 2.2 Continuous Variables

\begin{align}
X_{ijpt} &\geq 0 && \text{Flow of product } p \in \mathcal{P} \text{ from site } i \in \mathcal{I}  \to \text{ storage } j \in \mathcal{J} \text{ in period } t \in \mathcal{T^{harv}} \text{ [tonnes]} \\
X_{jkpp't} &\geq 0 && \text{Conversion flow of product } p \in \mathcal{P}\text{ to product } p'\in \mathcal{P}: p'\geq p \text{ from storage } j \to \text{ consumer } k \text{ in period } t  \in \mathcal{T^{harv}} \text{ [tonnes]}   \\
Y_{ipt} &\geq 0 && \text{Maximum harvest quantity of product } p \text{ at site } i  \text{ in period } t  \in \mathcal{T^{harv}} \text{ [tonnes]} \\
S_{jpt} &\geq 0 && \text{Inventory at storage } j \text{ in period } t  \in \mathcal{T}\text{ [tonnes]}
\end{align}

---

# 3. Parameters

| Parameter | Unit | Description | Typical Value |
|-----------|------|-------------|---------------|
| $\text{A}^{min}$ | yr | minimum  age | 3–5 years |
| $\text{A}^{max}$ | yr | maximum  age | 8–10 years |
| $\text{AREA}_i$ | ha | Site area | 8–20 ha |
| $c^{est}_i$ | €/ha | Establishment cost in € per ha| 10 €/ha |
| $c^{opp}_i$ | €/ha | Opportunity cost per ha | 5 – 15 €/ha |
| $c^{harv}_i$ | €/ha | Harvest and refitting cost per ha | 5 – 15 €/ha |
| $\eta_{pa}$ | t/ha | Base yield for product $p$, age $a$ | 4–60 t/ha |
| $R_{k,p}$ | €/t | Revenue for product $p$ by consumer $k$ | 80–250 € |
| $c^{tr-raw/pre}$ | €/t/km | Transport cost rates per ton and kilometer for raw and pre-treated biomass | 0.10/0.15 €/km-ton |
| $d_{ij}, d_{jk}$ | km | Distances | 30–200 km |
| $\text{CAP}^{stor}_j$ | t | Storage capacity | 300–500 t |
| $\text{CAP}^{proc}_j$ | t/yr | Processing capacity | 80–150 t/yr |
| $D^{max}_{kp}$ | t/yr | Demand at consumer $k$ | Varies |

---

# 4. Objective Function

\begin{align}
\max \, Z = \quad
&\sum_{k,p,t} R_{k,p} \cdot \sum_{j \in \mathcal{J}, p'\in \mathcal{P}: p' \leq p } X_{jkp'pt} && \text{Revenue} \\
&- \sum_{i} C^{est}_i \cdot \text{AREA}_i \cdot \sum_{t \in T} z_{i0t} && \text{Establishment} \\
&- \sum_{i} C^{opp}_i \cdot \text{AREA}_i \cdot \sum_{t \in T} (T^{max}-t) \cdot z_{i0t} && \text{Opportunity} \\
&- \sum_{i,(s,t)\in \mathcal{S}| s>0} C^{harv}_i \cdot \text{AREA}_i \cdot z_{ist} && \text{Harvest} \\
&- c^{tr-raw}\cdot\sum_{i,j,p,t} d_{ij} \cdot X_{ijpt} && \text{Transport (site→hub)} \\
&- c^{tr-pre}\cdot \sum_{j,k,p,t} d_{jk} \cdot X_{jkpt} && \text{Transport (hub→consumer)} \\
&- \sum_{j,p,t} c^{stor}_j \cdot S_{jpt} && \text{Storage}
\end{align}

---

# 5. Constraints 

## Constraint Set 1: AFS establishment

\begin{align}
\sum_{t \in \mathcal{T^+}} z_{i0t} &\leq 1 & \forall i \in \mathcal{I} \tag{C1}
\end{align}

**Interpretation:** If AFS is established in $t$, $z_{i,0,t}=1$.

---

## Constraint Set 2: Path connectivity

\begin{align}
\sum_{s = 0| (s,t) \in \mathcal{S}}^{t} z_{ist} &= \sum_{u = t+1| (s,t) \in \mathcal{S}}^{T_{max}+1} z_{itu} & \forall i \in \mathcal{I}, t \in \mathcal{T} \tag{C2}
\end{align}

**Interpretation:** If there is a harvest in $t$, there must be an predecessor in $s<t$ and successor in $u>t$.

---

## Constraint Set 3: Biomass yield calculation 

\begin{align}
Y_{ipt} &= \sum_{s=1}^{t-1} \eta_{p(t-s)} \cdot AREA_i \cdot z_{ist}  & \forall i \in \mathcal{I}, p \in \mathcal{P}, t \in \mathcal{T^{harv}} \tag{C3} \\
\end{align}

**Interpretation:** Biomass yield is restricted by age of plantage 

---


## Constraint Set 4: Yield Shipping  Linkage

\begin{align}
\sum_{j \in \mathcal{J}} X_{ijpt} &\leq Y_{ipt} & \forall i \in \mathcal{I}, p \in \mathcal{P}, t \in \mathcal{T^{harv}} \tag{C4}
\end{align}

**Interpretation:** Biomass shipped from site $i$ (product $p$, period $t$) is limited by available harvest.

---

## Constraint Set 5: Inventory Balance

\begin{align}
S_{jpt} &= S_{jp(t-1)} + \sum_{i} X_{ijpt} - \sum_{k, p' \geq p} X_{jkpp't} & \forall j \in \mathcal{J}, p \in \mathcal{P}, t \in \mathcal{T^{harv}} \tag{C5}
\end{align}

with $S_{j,p,0} = 0$.

**Interpretation:** Flow conservation at storage hubs.

---

## Constraint Set 6: Storage Capacity

\begin{align}
\sum_{p \in \mathcal{P}} S_{jpt} &\leq \text{CAP}^{stor}_j & \forall j \in\mathcal{J}, t \in \mathcal{T^{harv}} \tag{C6}
\end{align}

**Interpretation:** Total inventory limited by physical storage space.

---

## Constraint Set 7: Processing Capacity

\begin{align}
\sum_{i \in \mathcal{I}, p \in \mathcal{P}} X_{ijpt} &\leq \text{CAP}^{proc}_j & \forall j \in \mathcal{J}, t \in \mathcal{T^{harv}} \tag{C7}
\end{align}

**Interpretation:** Incoming material throughput limited by processing equipment.

---

## Constraint Set 8: Demand Satisfaction with Product Cascade

\begin{align}
\sum_{j \in \mathcal{J}, p'\in \mathcal{P}: p' \leq p } X_{jkp'pt} & \leq D^{max}_{kpt} & \forall k \in  \mathcal{K}, p \in  \mathcal{P}, t \in \mathcal{T^{harv}} \tag{C8}
\end{align}

**Interpretation:** Higher-quality products (lower $p$) can substitute for lower-quality demand.

**Cascade hierarchy:**
- Product 1 (chemical) → Satisfies any demand
- Product 2 (pulp) → Satisfies product 2 or 3 demand
- Product 3 (energy) → Only satisfies product 3 demand

---

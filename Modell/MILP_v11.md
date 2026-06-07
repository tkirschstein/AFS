---
title: "MILP Formulation for Agroforestry Biomass Supply Chain (v11)"
author: "SmartAgroforst"
date: "January 2026"
output:
  pdf_document: default
  html_document:
    df_print: paged
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
subtitle: "Multi-Cycle Harvesting with Age Lag Constraints"

header-includes:
  - \usepackage{booktabs}
  - \usepackage{longtable}
  - \usepackage{array}
  - \usepackage{multirow}
  - \usepackage{hyperref}
  - \usepackage{wrapfig}
  - \usepackage{colortbl}
  - \usepackage{pdflscape}
  - \usepackage{float}
  - \usepackage{xcolor}
  - \newcommand{\Tc}{\mathcal{T}}
  - \newcommand{\Pc}{\mathcal{P}}
  - \newcommand{\Qc}{\mathcal{Q}}
  - \newcommand{\Sc}{\mathcal{S}}
  - \newcommand{\Ic}{\mathcal{I}}
  - \newcommand{\Jc}{\mathcal{J}}
  - \newcommand{\Kc}{\mathcal{K}}
  - \newcommand{\problem}{AFS-SCD }
---


# Index Sets

\begin{align}
i \in \mathcal{I} &= \{1, \ldots, n_s\} && \text{Agroforestry sites} \\
j \in \mathcal{J} &= \{1, \ldots, n_j\} && \text{Processing facilities} \\
k \in \mathcal{K} &= \{1, \ldots, n_k\} && \text{Consumer sites} \\
t \in \mathcal{T} &= \{1, \ldots, T_{\max}\} && \text{Time periods} \\
t \in \mathcal{T}^{harv} &= \{A^{min}+1, \ldots, T_{\max}\} && \text{ Harvesting periods} \\
p \in \mathcal{P} &= \{1, 2, 3\} && \text{Product types} \\
(p,p') \in \mathcal{Q} &= \{(1,1),(1,2),(1,3),(2,2),(2,3),(3,3)\} && \text{Product conversion} \\
\Sc^{{est}}  &= \{(0,t) \mid t \in \{1,\dots,T^{\max}-A^{\min}\}\} &&\\
  \Sc^{{harv}} &= \{(s,t) \mid s,t \in \Tc \wedge t-s \in \{A^{\min},\dots,A^{\max}\}\}&&\\
  \Sc^{{term}} &= \{(s,T^{\max}+1) \mid s \in \{A^{\min}+1,\dots,T_{\max}\}\}&&\\
  \Sc               &= \Sc^{{est}} \cup \Sc^{{harv}} \cup \Sc^{{term}}&&
\end{align}

---

# Decision Variables

\begin{align}
z_{ist} & \in [0, AREA_i] && \text{size of AFS of site} i \text{used in periods} t,s \text{ with } (s,t) \in \mathcal{S}\\
X_{ijpt} &\geq 0 && \text{Flow of product } p \in \mathcal{P} \text{ from } i   \to  j  \text{ in } t \in\mathcal{T}^{harv}  \\
X_{jkpp't} &\geq 0 && \text{Conversion flow of } p \text{ to  } p' \text{ from } j \to k \text{ in } t  \in \mathcal{T}^{harv}    \\
S_{jpt} &\geq 0 && \text{Inventory at storage } j \text{ in period } t  \in \mathcal{T}
\end{align}

---

# Objective Function

\begin{align}
\max\; Z \;=\;&
  \sum_{k \in \Kc}\sum_{p \in \Pc} R_{kp} \cdot \sum_{j \in \Jc, }\sum_{(p',p) \in \Qc} \sum_{t \in \Tc}  X_{jkp'pt}  \label{eq:obj-revenue} \\
&- \sum_{i \in \Ic} c^{{est}}_i \cdot  \sum_{(0,t) \in \Sc^{{est}}} z_{i0t} \notag \\
&- \sum_{i \in \Ic} (c^{main}_{i} + c^{{opp}}_i) \cdot (t-s) \cdot \sum_{(s,t) \in \Sc^{{harv}}} z_{ist} \notag \\
&- \sum_{i \in \Ic} c^{{harv}}_i  \cdot \sum_{(s,t) \in \Sc^{{harv}}} z_{ist} \notag \\
&- \sum_{i \in \Ic}\sum_{j \in \Jc}\sum_{p \in \Pc}\sum_{t \in \Tc} c^{\text{tr-raw}}_p \cdot d_{ij} \cdot X_{ijpt} \notag \\
&- \sum_{j \in \Jc}\sum_{k \in \Kc}\sum_{(p,p') \in \Qc} \sum_{t \in \Tc} c^{\text{tr-pre}}_{p'} \cdot d_{jk} \cdot X_{jkpp't} \notag \\
&- \sum_{j \in \Jc}\sum_{p \in \Pc}\sum_{t \in \Tc}  c^{{stor}}_j \cdot S_{jpt}  \notag
\end{align}

---

# Constraints 

\begin{align}
  % --- C1: Establishment (at most once per site) ---
  \sum_{(0,t) \in \Sc^{{est}}} z_{i0t} &\leq AREA_i && \forall\; i \in \Ic \label{eq:C1} \\[2pt]
  %
  % --- C2: Path connectivity (flow conservation) ---
  \sum_{(s,t) \in \Sc} z_{ist} &= \sum_{(t,u) \in \Sc} z_{itu} && \forall\; i \in \Ic,\; t \in \Tc \label{eq:C2} \\[2pt]
  %
  % --- C3: Biomass yield from Gompertz coefficients ---
  \sum_{j \in \Jc} X_{ijpt}    &\leq \sum_{(s,t) \in \Sc^{{harv}}} \eta_{p(t-s)} \cdot z_{ist} && \forall\; i \in \Ic,\; p \in \Pc,\; t \in \Tc^{{harv}}     \label{eq:C3} \\[2pt]
  %
  % --- C5: Inventory balance ---
  S_{jpt}    = S_{jp,t-1} + &\sum_{i \in \Ic} X_{ijpt} - \sum_{k \in \Kc,\,(p,p') \in \Qc} X_{jkpp't}   && \forall\; j \in \Jc,\; p \in \Pc,\; t \in \Tc^{harv}  \label{eq:C5} \\[2pt]
  %
  % --- C6: Storage capacity ---
  \sum_{p \in \Pc} S_{jpt}    &\leq \text{CAP}^{{stor}}_j    && \forall\; j \in \Jc,\; t \in \Tc     \label{eq:C6} \\[2pt]
  %
  % --- C7: Processing capacity ---
  \sum_{i \in \Ic,\, p \in \Pc} X_{ijpt}  &\leq \text{CAP}^{{proc}}_j && \forall\; j \in \Jc,\; t \in \Tc^{harv} \label{eq:C7} \\[2pt]
  %
  % --- C8: Demand satisfaction with quality cascade ---
  \sum_{j \in \Jc,\,(p',p) \in \Qc} X_{jkp'pt} &\leq D^{\max}_{kpt} && \forall\; k \in \Kc,\; p \in \Pc,\; t \in \Tc^{harv}     \label{eq:C8} \\[2pt]
  %
  % --- Domain: continuous non-negativity ---
  z_{ist},\; X_{ijpt},\; &X_{jkp'pt},\; S_{jpt} \geq 0 && \forall\; i,j,k,p,p',t     \label{eq:D2}
  %
\end{align}

---

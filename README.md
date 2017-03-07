# PlogP_TRN
Reconstruct a Transcriptional Regulatory Network using the principle of Maximum Entropy.

## About

PlogP_TRN implements a regulatory network reconstruction method based on a method proposed by [T.R. Lezon et al. Using the principle of entropy maximization to infer genetic interaction networks from gene expression patterns, PNAS, 2006, 103(50):19033-19038 doi:10.1073/pnas.0609152103](http://doi.org/10.1073/pnas.0609152103).

Briefly, the method is based on finding <a href="https://www.codecogs.com/eqnedit.php?latex=\rho&space;(\vec{x})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\rho&space;(\vec{x})" title="\rho (\vec{x})" /></a> where <a href="https://www.codecogs.com/eqnedit.php?latex=\vec{x}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\vec{x}" title="\vec{x}" /></a> are the observed gene expression values by maximizing the Shannon entropy:

<a href="https://www.codecogs.com/eqnedit.php?latex=S=-\sum_{\vec{x}}\rho&space;(\vec{x})ln\rho&space;(\vec{x})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?S=-\sum_{\vec{x}}\rho&space;(\vec{x})ln\rho&space;(\vec{x})" title="S=-\sum_{\vec{x}}\rho (\vec{x})ln\rho (\vec{x})" /></a> 

which leads to a Boltzmann(-ish) distribution:

<a href="https://www.codecogs.com/eqnedit.php?latex=\rho&space;(\vec{x})\propto&space;e^{-\frac{\sum_{ij}x_{i}M_{ij}x_{j}}{2}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\rho&space;(\vec{x})\propto&space;e^{-\frac{\sum_{ij}x_{i}M_{ij}x_{j}}{2}}" title="\rho (\vec{x})\propto e^{-\frac{\sum_{ij}x_{i}M_{ij}x_{j}}{2}}" /></a>

where <a href="https://www.codecogs.com/eqnedit.php?latex=M_{ij}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?M_{ij}" title="M_{ij}" /></a> is the "interaction strength" between genes *i* and *j*.

## Installation

Download the 

Requirements:

* [Julia](http://julialang.org)
* [IJulia](https://github.com/JuliaLang/IJulia.jl)

## Usage



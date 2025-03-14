## README

This repository contains R codes implementing methods developed in my thesis for curve fitting and derivative estimation in longitudinal data using P-splines.

### Repository Structure

### SPNLM Models

This folder includes:

* ***Implementation of SPNLME***: Scripts for fitting Semiparametric nonlinear mixed-effects (SPNLME) models. The Preece-Baines model is used as parametric component and  P-splines model as nonparamteric component.
* ***Classification Methods***: Scripts for applying Discriminant classification  using SPNLME

###  SITAR Modifications

This folder includes implementations of SITAR model(Cole et al., 2010) modifications:

* ***SITAR with P-Splines***: Modification of the SITAR model using penalized splines.
* ***SITAR with SAEM Algorithm***: Implementation of the Stochastic Approximation Expectation-Maximization (SAEM) algorithm for parameter estimation in SITAR model.
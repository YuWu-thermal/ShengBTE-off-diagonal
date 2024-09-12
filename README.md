# An extension module to ShengBTE for computing glass-like thermal conductivity

## $\color{red}{\rm{Attention}}$：

An issue in the phonon_routines.f90 file uploaded earlier have been updated on date 2024/6/24
Additional output can be obtained: BTE.cumulative_Glass_kappaVsOmega_tensor, BTE.kappa_glass_resolved 2024/9/12

## Based version:

FourPhonon1.0

## Output:

1. BTE.Kappa_Glass_like_RTA
2. Cumulative glass-like thermal conductivity: BTE.cumulative_Glass_kappaVsOmega_tensor  
3. Resolved glass-like thermal conductivity associated with different phonon pairs: BTE.kappa_glass_resolved   Format:(omega1, omega2, kappaxx, kappaxy...) Visualization: python plot_resolved_kappa.py
   
## Authors and references:

Yu Wu: wuyu9573@qq.com

Reference:

Yu Wu et al.,  [J. Phys. Chem. Lett. 2023, 14, 11465−11473](https://doi.org/10.1021/acs.jpclett.3c02940)

Yu Wu et al.,  [Phys. Rev. B 2024, 109, 214307](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.109.214307)

## The original ShengBTE:

https://github.com/FourPhonon/FourPhonon


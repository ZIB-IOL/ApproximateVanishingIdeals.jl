# Border construction
This section contains information about the construction of borders for $\texttt{OAVI}$.

## Initial construction
The first time constructing the border we need to handle non-unique terms, purge them and then only continue with the relevant terms and corresponding evaluations.
```@docs
AVI.construct_border
```

```@docs
AVI.purge
```

## Reconstruction of the border
When applying the transformation $\mathcal{G}$ found by OAVI, we need to reconstruct the border. This can be sped up by using the computations done in `construct_border` to avoid unnecessary recomputations of known values. 
```@docs
AVI.reconstruct_border
```

## Index
```@index
Pages=["3_border.md"]
```

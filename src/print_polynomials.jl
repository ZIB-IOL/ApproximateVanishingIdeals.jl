"""
Prints the polynomials obtained through OAVI as a LaTeX string.
"""
function print_polynomials(sets; digits::Int64=2, render=false)
    deg_ids = sets.O_degree_indices

    lt_idx_start = 1
    lt_idx_end = 0

    all_polys = []

    for i in eachindex(sets.G_coefficient_vectors)
        coeff_vec = sets.G_coefficient_vectors[i]
        if coeff_vec !== nothing
            lt_idx_end += size(coeff_vec, 2)

            coeff_vec = coeff_vec[end:-1:1, :]

            # shift -1 due to O_deg structure
            if i <= length(sets.O_degree_indices)
                idx = deg_ids[i]
            else
                idx = size(sets.O_terms, 2)
            end

            # get only non-leading terms necessary for current coefficient vector
            terms = sets.O_terms[:, 1:idx]

            # get leading terms of current degree
            lts = sets.leading_terms[:, lt_idx_start:lt_idx_end]
            
            # all terms together in deg-lex ordering to fit the coefficient vector
            all_terms = hcat(lts, terms)
            all_terms, _ = AVI.deg_lex_sort(all_terms)
            all_terms = all_terms[:, end:-1:1]
    
            for (j, col) in enumerate(eachcol(coeff_vec))   
                # build polynomial 
                poly_string = ""
                for (k, coeff) in enumerate(col)
                    coeff = round(coeff, digits=digits)
                    if coeff != 0.0
                        # omit coefficient if coeff = 1
                        if coeff == 1.0
                            # bring term into LaTeX format
                            poly_string *= convert_term_to_latex(all_terms[:, k])
                            poly_string *= " + "

                        # add coeff and term if > 0
                        elseif coeff > 0.0
                            poly_string *= string(coeff)
                            poly_string *= convert_term_to_latex(all_terms[:, k])
                            poly_string *= " + "

                        # replace '+' by '-' if negative
                        elseif coeff < 0.0
                            poly_string = poly_string[1:length(poly_string)-2]
                            poly_string *= "- "
                            poly_string *= string(abs(coeff))
                            poly_string *= convert_term_to_latex(all_terms[:, k])
                            poly_string *= " + "
                        end
                    else
                        continue
                    end
                end

                # append only if string not empty
                if poly_string !== ""
                    poly_string = poly_string[1:length(poly_string)-3]
                    append!(all_polys, [poly_string])
                end
            end
            lt_idx_start = lt_idx_end + 1
        end 
    end
    for (i, poly) in enumerate(all_polys)
        println(poly)
    end
end


"""Converts a term given as a vector denoting the exponents into LaTeX string"""
function convert_term_to_latex(term::Vector{Int64})
    term_string = ""
    for i in eachindex(term)
        if term[i] == 1
            term_string *= "x_$(i)"
        elseif term[i] > 0
            term_string *= "x_{$(i)}^{$(term[i])}"
        end
    end
    if term_string == ""
    end
    return term_string
end

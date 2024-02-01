# Transforming data using $\texttt{VCA}$
If you whis to use $\texttt{VCA}$ instead of $\texttt{OAVI}$ to obtain your transformation, then here we will briefly explain how you can go about it.

## Choosing the right transformation
Unlike $\texttt{OAVI}$, $\texttt{VCA}$ does not enjoy a known relation between the vanishing parameter $\psi$ and the degree of constructed polynomials. This makes a wider range of $\psi$'s to compare necessary, even if one has some knowledge about which polynomial degrees may prove better. Luckily, the way to go about finding the best transformation is very similar to the one presented for $\texttt{OAVI}$.

Once again, let's begin by creating or reading in our data set.

````julia docs_vca_transformation
using AVI
using DataFrames
using CSV
using Random

# create data
X = rand(100, 3)

# put it into a DataFrame
df = DataFrame(X, :auto)

# write it to 'example.csv' in 'examples' folder
CSV.write("examples/example_vca.csv", df);

# read the data
read_df = CSV.read("examples/example_vca.csv", DataFrame, types=Float64)

# convert DataFrame to Matrix
data = Matrix(read_df)
````

Just as with $\texttt{OAVI}$, `fit_vca` expects the data matrix to be of type `Float64`.

As with the previous example, the next step is to choose some range of $\psi$'s you wish to test. We recommend to choose a wider range, for example `psis=[0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001]`. Now, just as with $\texttt{OAVI}$, we wish to test and compare the transformations for different values of $\psi$.

````julia docs_vca_transformation
# split train and test data
X_train, X_test = data[1:90, :], data[91:end, :]

# define psis
psis=[0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001]

# init best transform and best sets
best_transform_vca, best_sets_vca = nothing, nothing

# init error
error = Inf64

# loop psis
for psi in psis
    # transformation with current psi
    X_train_transformed_vca, sets_train_vca = fit_vca(X_train; psi=psi)

    """
    if X_train_transformed_vca !== nothing
        new_error = classifier(X_train_transformed_vca, class_labels)
        if new_error < error
            error = new_error
            best_transform_vca, best_sets_vca = X_train_transformed_vca, sets_train_vca
        end
    end
    """
end
````

Use your favorite classification method to compare the transformations and choose the one that fits your data and problem best.

## Applying the transformation
After having found the transformation best suited to your data, you can apply this transformation to the data you wish to classify. Passing `best_sets_vca` and `X_test` to `apply_V_transformation` applies the $\texttt{VCA}$ transformation stored in `best_sets_vca` to `X_test`. We compute some `best_transform_vca` and 'best_sets_vca' as to not get an error.

````julia docs_vca_transformation
# best transform and best sets for VCA
best_transform_vca, best_sets_vca = fit_vca(X_train; psi=0.001)

# apply VCA transformation to test set
X_test_transformed_vca, sets_test_vca = apply_V_transformation(best_sets_vca, X_test)
````

Once again, the data to which the transformation is applied should have the same second dimension as `X_train`, which means `size(X_train, 2) == size(X_test, 2)` should hold.

We perform the same sanity check as for $\texttt{OAVI}$ to see if the dimension of the transformation matches.

````julia docs_vca_transformation
# get all polynomials into single matrix
polys = AVI.V_to_matrix(best_sets_vca)

println("Number of vanishing polynomials: ", size(polys, 2))
size(best_transform_vca, 2) == size(X_test_transformed_vca, 2) == size(polys, 2)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*


array[] real divide_vector(vector vec, array[] real c) {
    int n = num_elements(vec);
    array[n] real product;
    for (i in 1:n) {
        product[i] = vec[i] / c[i];
    }
    return product;
}

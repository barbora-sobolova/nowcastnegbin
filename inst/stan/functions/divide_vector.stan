array[] real divide_vector(vector vec, array[] real c) {
    int n = num_elements(vec);
    if (num_elements(c) != n) {
        reject("divide_vector: `vec` and `c` must have the same length.");
    }
    array[n] real product;
    for (i in 1:n) {
      if (c[i] <= 0) {
        reject("divide_vector: divisor must be > 0 at index ", i);
      }
      product[i] = vec[i] / c[i];
    }
    return product;
}

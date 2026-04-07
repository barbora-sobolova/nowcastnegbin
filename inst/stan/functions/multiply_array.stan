array[] real multiply_array(real c, array[] real arr) {
    int n = num_elements(arr);
    array[n] real product;
    for (i in 1:n) {
        product[i] = c * arr[i];
    }
    return product;
}

array[] real multiply_array(array[] real c, array[] real arr) {
    int n = num_elements(arr);
    if (num_elements(c) != n) {
        reject("multiply_array: `c` and `arr` must have the same length.");
    }
    array[n] real product;
    for (i in 1:n) {
        product[i] = c[i] * arr[i];
    }
    return product;
}

array[] real multiply_array(vector c, array[] real arr) {
    int n = num_elements(arr);
    if (num_elements(c) != n) {
        reject("multiply_array: `c` and `arr` must have the same length.");
    }
    array[n] real product;
    for (i in 1:n) {
        product[i] = c[i] * arr[i];
    }
    return product;
}

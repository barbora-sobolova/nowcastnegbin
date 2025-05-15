array[] real multiply_array_by_vector(vector vect, array[] real arr) {
    int n = num_elements(arr);
    array[n] real product;
    for (i in 1:n) {
        product[i] = vect[i] * arr[i];
    }
    return product;
}

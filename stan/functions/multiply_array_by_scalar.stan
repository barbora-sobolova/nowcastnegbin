array[] real multiply_array_by_scalar(real scal, array[] real arr) {
    int n = num_elements(arr);
    array[n] real product;
    for (i in 1:n) {
        product[i] = scal * arr[i];
    }
    return product;
}

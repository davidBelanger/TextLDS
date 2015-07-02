function assertEqual(A,B)
    assert(all(abs(A(:) - B(:)) < 0.0000001));
    
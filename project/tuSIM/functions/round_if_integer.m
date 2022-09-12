function x = round_if_integer(x, err_msg)
    assert(abs(round(x)-x) < 1e-6, err_msg)
    x = round(x);
end
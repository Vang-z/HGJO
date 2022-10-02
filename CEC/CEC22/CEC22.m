function [lb, ub, dim, f_obj, f_name] = CEC22(i_func)
    switch i_func
        case 1
            f_obj = @F1;
            lb = -100;
            ub = 100;
            dim = 10;
            f_name = "Zakharov Function";
        case 2
            f_obj = @F2;
            lb = -100;
            ub = 100;
            dim = 10;
            f_name = "Rosenbrock's Function";
        case 3
            f_obj = @F3;
            lb = -100;
            ub = 100;
            dim = 10;
            f_name = "Schaffer's F7";
        case 4
            f_obj = @F4;
            lb = -100;
            ub = 100;
            dim = 10;
            f_name = "Rastrigin's Function";
        case 5
            f_obj = @F5;
            lb = -100;
            ub = 100;
            dim = 10;
            f_name = "Levy Function";
        case 6
            f_obj = @F6;
            lb = -100;
            ub = 100;
            dim = 10;
            f_name = "Hybrid Function 1";
        case 7
            f_obj = @F7;
            lb = -100;
            ub = 100;
            dim = 10;
            f_name = "Hybrid Function 2";
        case 8
            f_obj = @F8;
            lb = -100;
            ub = 100;
            dim = 10;
            f_name = "Hybrid Function 3";
        case 9
            f_obj = @F9;
            lb = -100;
            ub = 100;
            dim = 10;
            f_name = "Composition Function 1";
        case 10
            f_obj = @F10;
            lb = -100;
            ub = 100;
            dim = 10;
            f_name = "Composition Function 2";
        case 11
            f_obj = @F11;
            lb = -100;
            ub = 100;
            dim = 10;
            f_name = "Composition Function 3";
        case 12
            f_obj = @F12;
            lb = -100;
            ub = 100;
            dim = 10;
            f_name = "Composition Function 4";
    end
end


function o = F1(x)
    o = cec22_test_func(x', 1);
end

function o = F2(x)
    o = cec22_test_func(x', 2);
end

function o = F3(x)
    o = cec22_test_func(x', 3);
end

function o = F4(x)
    o = cec22_test_func(x', 4);
end

function o = F5(x)
    o = cec22_test_func(x', 5);
end

function o = F6(x)
    o = cec22_test_func(x', 6);
end

function o = F7(x)
    o = cec22_test_func(x', 7);
end

function o = F8(x)
    o = cec22_test_func(x', 8);
end

function o = F9(x)
    o = cec22_test_func(x', 9);
end

function o = F10(x)
    o = cec22_test_func(x', 10);
end

function o = F11(x)
    o = cec22_test_func(x', 11);
end

function o = F12(x)
    o = cec22_test_func(x', 12);
end

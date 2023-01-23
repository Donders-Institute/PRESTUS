classdef Tests < matlab.unittest.TestCase

    methods(Test)
        function load(testCase)
            tests = {
                % YAML | expected result
                "test # comment", "test"
                "1.23", 1.23
                "True", true
                "1", 1
                "[1, 2]", {1, 2}
                "[1, 2, True, test]", {1, 2, true, "test"}
                "{}", struct()
                sprintf("12!: 1\n12$: 2"), struct("x12_", 1, "x12__1", 2)
                sprintf("a: test\nb: 123\nc:\n  d: test2\n  e: False"), struct("a", "test", "b", 123, "c", struct("d", "test2", "e", false))
                ".nan", NaN
                ".inf", inf
                "-.inf", -inf
                "null", yaml.Null
                "", yaml.Null
                "~", yaml.Null
                "2019-09-07T15:50:00", datetime(2019, 9, 7, 15, 50, 0, "TimeZone", "UTC")
                "2019-09-07 15:50:00", datetime(2019, 9, 7, 15, 50, 0, "TimeZone", "UTC")
                "2019-09-07", datetime(2019, 9, 7, "TimeZone", "UTC")
                "2019 09 07 15:50:00", "2019 09 07 15:50:00"
                };

            for test = tests'
                [s, expected] = test{:};
                actual = yaml.load(s);
                testCase.verifyEqual(actual, expected);
            end
        end

        function load_converToArray(testCase)
            tests = {
                % YAML | expected result
                "[1]", 1
                "[1, 2]", [1, 2]
                "[[1, 2], [3, 4]]", [1, 2; 3, 4]
                "[[[1, 2], [3, 4]], [[5, 6], [7, 8]]]", {[1, 2; 3, 4], [5, 6; 7, 8]}
                "[[1, 2], [3]]", {[1, 2], 3}
                "[]", []
                "[[1, 2], []]", {[1, 2], []}
                "[1, true]", {1, true}
                "[[a, b], [c, d]]", ["a", "b"; "c", "d"]
                "[null, 1]", {yaml.Null, 1}
                "[null, null]", [yaml.Null, yaml.Null]
                };

            for test = tests'
                [s, expected] = test{:};
                actual = yaml.load(s, "ConvertToArray", true);
                testCase.verifyEqual(actual, expected);
            end
        end

        function dump(testCase)
            tests = {
                % Data | expected YAML
                "test", "test"
                'test', "test"
                't', "t"
                1.23, "1.23"
                int32(1), "1"
                true, "true"
                struct("a", "test", "b", 123), "{a: test, b: 123.0}"
                {}, "[]"
                {1, "test"}, "[1.0, test]"
                {1; "test"}, "[1.0, test]"
                {1, {2, 3}}, sprintf("- 1.0\n- [2.0, 3.0]")
                nan, ".NaN"
                inf, ".inf"
                -inf, "-.inf"
                yaml.Null, "null"
                [1, 2], "[1.0, 2.0]"
                ["a", "b"], "[a, b]"
                [true, false], "[true, false]"
                [], "[]"
                zeros(1, 0), "[]"
                zeros(0, 1), "[]"
                int32(ones(2, 2, 2)), sprintf("- - [1, 1]\n  - [1, 1]\n- - [1, 1]\n  - [1, 1]")
                num2cell(int32(ones(2, 2, 2))), sprintf("- - [1, 1]\n  - [1, 1]\n- - [1, 1]\n  - [1, 1]")
                int32(ones(2, 1, 2)), sprintf("- [1, 1]\n- [1, 1]")
                };

            for test = tests'
                [data, expected] = test{:};
                expected = expected + newline;
                actual = yaml.dump(data);
                testCase.verifyEqual(actual, expected);
            end
        end

        function dump_3dcell(testCase)
            data = num2cell(ones(2, 2, 2));
            data{1, 1, 2} = "a";
            data{1, 2, 1} = yaml.Null;
            data{2, 1, 1} = {1, 2};

            expected = sprintf("- - [1.0, a]\n  - ['null', 1.0]\n- - - [1.0, 2.0]\n    - 1.0\n  - [1.0, 1.0]\n");

            actual = yaml.dump(data);
            testCase.verifyEqual(actual, expected);
        end

        function dump_unsupportedTypes(testCase)
            tests = {
                % Data | expected error
                ones(2, 2, 2, 2), "yaml:dump:HigherDimensionsNotSupported"
                num2cell(ones(2, 2, 2, 2)), "yaml:dump:HigherDimensionsNotSupported"
                datetime(2022, 2, 13), "yaml:dump:TypeNotSupported"
                "test $%&? adfasdf", "yaml:dump:NullPlaceholderNotAllowed"
                };

            for test = tests'
                [data, errorId] = test{:};
                func = @() yaml.dump(data);
                testCase.verifyError(func, errorId);
            end
        end

        function dump_style(testCase)
            data.a = 1;
            data.b = {3, {4}};

            tests = {
                "block", sprintf("a: 1.0\nb:\n- 3.0\n- - 4.0\n")
                "flow", sprintf("{a: 1.0, b: [3.0, [4.0]]}\n")
                "auto", sprintf("a: 1.0\nb:\n- 3.0\n- [4.0]\n")
                };

            for iTest = 1:size(tests, 1)
                [style, expected] = tests{iTest, :};
                actual = yaml.dump(data, style);
                testCase.verifyEqual(actual, expected);
            end

        end

        function dumpFile(testCase)
            data = struct("a", 1.23, "b", "test");
            expected = "{a: 1.23, b: test}";
            if ispc
                expected = expected + sprintf("\r\n");
            else
                expected = expected + sprintf("\n");
            end

            testPath = tempname;

            yaml.dumpFile(testPath, data)
            fid = fopen(testPath);
            actual = string(fscanf(fid, "%c"));
            fclose(fid);

            testCase.verifyEqual(actual, expected);

            delete(testPath)
        end

        function loadFile(testCase)
            data = struct("a", 1.23, "b", "test");

            testPath = tempname;
            yaml.dumpFile(testPath, data)
            actual = yaml.loadFile(testPath);

            testCase.verifyEqual(actual, data);
            delete(testPath)
        end

        function loadFile_convertToArray(testCase)
            data = {1, 2};
            expected = [1, 2];

            testPath = tempname;
            yaml.dumpFile(testPath, data)
            actual = yaml.loadFile(testPath, "ConvertToArray", true);

            testCase.verifyEqual(actual, expected);
            delete(testPath)
        end

        function isNull(testCase)

            testCase.verifyTrue(yaml.isNull(yaml.Null))

            nonNulls = {NaN, missing, "", "a", datetime(2022, 1, 1), NaT, '', {}, [], inf, -inf, 1};

            for i = 1:length(nonNulls)
                testCase.verifyFalse(yaml.isNull(nonNulls{i}))
            end
        end

        function null(testCase)
            tests = {
                % Input arguments | Expected
                {}, yaml.Null
                {1}, yaml.Null
                {2}, repmat(yaml.Null, 2, 2)
                {2, 3}, repmat(yaml.Null, 2, 3)
                {3, 2}, repmat(yaml.Null, 3, 2)
                };

            for test = tests'
                [args, expected] = test{:};
                actual = yaml.Null(args{:});
                testCase.assertEqual(actual, expected)
            end
        end
    end
end

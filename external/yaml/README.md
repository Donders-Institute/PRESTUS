# yaml
[YAML](https://yaml.org/) parser and emitter for MATLAB R2019b or newer. Based on [SnakeYAML](https://bitbucket.org/snakeyaml/snakeyaml/src/master/).

[![View yaml on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/106765-yaml)

## Examples
### Load and dump
```Matlab
>> data.a = [1.23, 4.56];
>> data.b = {int32(2), {true, "hello", yaml.Null}};

>> s = yaml.dump(data)
    "a: [1.23, 4.56]
     b:
     - 2
     - [true, hello, null]
     "
   
>> result = yaml.load(s)
    a: {[1.2300]  [4.5600]}
    b: {[2]  {1×3 cell}}
```

### Read and write files
```Matlab
>> yaml.dumpFile("test.yaml", data)
>> result = yaml.loadFile("test.yaml")
    a: {[1.2300]  [4.5600]}
    b: {[2]  {1×3 cell}}
```

### Styles
```Matlab
>> s = yaml.dump(data, "auto")  % default
    "a: [1.23, 4.56]
     b:
     - 2
     - [true, hello, null]
     "
     
>> s = yaml.dump(data, "block")
    "a: 
     - 1.23
     - 4.56
     b:
     - 2
     - - true
       - hello
       - null
     "
     
>> s = yaml.dump(data, "flow")
    "{a: [1.23, 4.56], b: [2, [true, hello, 'null']], c: [2, [true, hola]]}
     "
```
### YAML null
```Matlab
>> result = yaml.load("null")
    Null
    
>> yaml.isNull(result)
   1
   
>> s = yaml.dump(yaml.Null)
    "null
     "
```

### Load YAML sequences as MATLAB standard arrays
By default, sequences are loaded as nested cell arrays to distinguish between YAML scalars and YAML one-element sequences and to supported mixed type sequences. If you use the `ConvertToArray` option, sequences are converted to 1D or 2D standard arrays if possible:
```Matlab
>> yaml.load("[[1, 2], [3, 4]]", "ConvertToArray", true)
     1     2
     3     4
```

### Control dumping behaviour for MATLAB arrays
Since every MATLAB scalar is always an array and every array technically has at least 2 dimensions, there exists two ambiguities when dumping arrays:
- *MATLAB scalar* &rarr; *YAML scalar* (default) or *YAML one-element sequence*
- *MATLAB vector* &rarr; *YAML sequence* (default) or *YAML sequence containing one YAML sequence*

To avoid theses ambiguities and get consistent conversion behaviour, convert all your array data to **nested vector cells** before dumping them.
```Matlab
>> yaml.dump({1})
    "[1.0]
    "
>> yaml.dump({{1, 2}})
    "- [1.0, 2.0]
    "
```


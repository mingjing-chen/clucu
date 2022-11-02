# 调试的方法

1. 创建静态库
```bash
gcc -c  clucu_error.c clucu_utils.c clucu_core.c clucu_background.c clucu_massfunction.c clucu_survey.c clucu_probe.c

ar –crv libclucu.a  clucu_error.o clucu_utils.o clucu_core.o clucu_background.o clucu_massfunction.o clucu_survey.o clucu_probe.o

mv libclucu.a ../lib
```
修改`tasks.json`文件
```json
"args":[
"-I/home/mjchen/Program/cluster/newcode/include",
"-L/home/mjchen/Program/cluster/newcode/lib",
"-lclucu",
]
```

2. 创建动态库
```bash
gcc –fpic  –shared –o libclucu.so  clucu_error.c clucu_utils.c clucu_core.c clucu_background.c clucu_massfunction.c clucu_survey.c clucu_probe.c
```
修改`tasks.json`文件
```json
"args":[
"-I/home/mjchen/Program/cluster/newcode/include",
"/home/mjchen/Program/cluster/newcode/source/libclucu.so",
]
```

3. 直接在`tasks.json`文件中加入源文件`clucu_*.c`路径
```json
"args":[
"-I/home/mjchen/Program/cluster/newcode/include",
"/home/mjchen/Program/cluster/newcode/source/clucu_*.c",
]
```

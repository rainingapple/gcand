# gcand

Source code for **gcand.so** metioned in our paper.

**gcand** simplify the LKH2 algorithm to generate candidatesets for GNN model. 

We use **Ctype** as the bridge between C and Python , use Function parameter to refrain from the expensive IO cost

## Changes

- Using Function parameter to get TSP Data instead of Reading from files
- Fixed EDGE_WEIGHT_TYPE to EUC_2D , COORDINATE_TYPE to NODE_COORD_SECTION
- Only retain the CandidateSet part (CreateCandidateSet.c) in LKH , remove others
- No main fuction , Only LKHmain for python to invoke

## Usage

### Download

```shell
git clone https://github.com/rainingapple/gcand.git
```

### Build

Using gcc -shared to build a .so file

```shell
cd gcand
gcc -I SRC/INCLUDE SRC/*.c -fPIC -shared -o gcand.so
```

### Load and Call

Move file to corresponding folder

```shell
mv gcand.so "your folder path"
```

Load Library and call LKHmain fuction as fellow

```python
cand_genrater = cdll.LoadLibrary('./gcand.so')
cand_genrater.LKHmain.restype = c_char_p
cand = cand_genrater.LKHmain(create_string_buffer(data_str.encode('utf-8'))).decode("utf-8").split("\n")
```

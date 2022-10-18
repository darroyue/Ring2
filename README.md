# Ring

Repository for the source code of the engine presented in the paper Worst-case Optimal Graph Joins in Almost No Space.

## Instructions

To run our code, **we have to install an extended version of the library SDSL**. You can find the library and the instructions to install it in our webpage.

After the extended version of SDSL is installed, we have to clone this repository and follow these steps:

1. Create our `build` folder and compile the code:
```Bash
mkdir build
cd build
cmake ..
make
```

Check that there is no errors.

2. Download the version of Wikidata that you want to use:

- [Wikidata Filtered (about 80M triples)](http://compact-leapfrog.tk/files/wikidata-filtered-enumerated.dat).
- [Wikidata (about 1000M triples)](http://compact-leapfrog.tk/files/wikidata-enumerated.dat.gz). Note that this file is compressed.

Now put the .dat file inside a folder.

3. Building the index. After compiling the code we should have an executable called `build-index` in `build`. Now run:

```Bash
./build-index <absolute-path-to-the-.dat-file> <type-ring>
```

`<type-ring>` can take two values: `ring` or `c-ring`. Both are implementations of our ring index but using plain and compressed bitvectors, respectively.
This will generate the index in the folder where the `.dat` file is located. The index is suffixed with `.ring` or `.c-ring` according to the second argument.

4. Querying the index. In `build` folder, you should find another executable file called `query-index`. To solve the queries you should run:

```Bash
./query-index <absoulute-path-to-the-index-file> <absolute-path-to-the-query-file>
```

Note that the second argument is the path to a file that contains all the queries. The queries of our benchmark are in `Queries`:

- The file `Queries-wikidata-benchmark.txt` can be run with `wikidata-filtered-enumerated.dat`.
- The file `Queries-bgps-limit1000.txt` contains the queries of `wikidata-enumerated.dat`.

After running that command, you should see the number of the query, the number of results, and the elapsed time of each one of the queries with the following format:
```Bash
<query number>;<number of results>;<elapsed time>
```
---

At the moment, we can find the rest of the complementary material at [this webpage](http://compact-leapfrog.tk/). Note that we will find instructions to run the code there, and although the instructions are different from the ones in this repository, they should work too.

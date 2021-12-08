# Assingment 3 - Carl Holmberg

I solved part 1.
Part 2 gave me some trouble. I tried many things but my code never seamed to work. I will not include the code at all.
Therefore the rest of this document is written from the perspective of only answering part 1.

## How to run

You need to have python and pip installed to run this program.

Install the required packages by running

```
$ pip install -r requirements.txt
```

Thereafter to run the program you just type

```
./identify_order.py [FILE]
```

where `[FILE]` is the input file.

## Approach and assumptions

The approach taken on part 1 is to start by finding the end nodes. This is done by finding the alpha-carbon atoms with only one neighbor. From that the chain is found by repeatatly looking for the next atom in the chain.

There was assumed to only be one chain.
An assumption was also made that there is no "overlap" of two atoms, for example that an atom has more than two neighbours.

## Output

_part 1 - test_q1.txt_

```
❯ ./identify_order.py data/test_q1.txt
2
4
1
3
5
6
7
9
10
8
Total alpha carbon atoms: 10
```

![output image from test_q1.txt](https://i.imgur.com/1zQED5P.png)
_part 1 - data_q1.txt_

```
❯ ./identify_order.py data/data_q1.txt
2
1
10
5
7
6
4
3
9
8
Total alpha carbon atoms: 10
```

![output image from data_q1.txt](https://i.imgur.com/NwVsIcq.png)

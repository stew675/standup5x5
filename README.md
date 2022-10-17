### Stand Up Maths 5 x 5 unique letters solutions

My solutions to the problem posed by Matt Parker of Stand Up Maths

See: https://www.youtube.com/watch?v=_-AfhLQfb6w


### Building and Running

Just run make and execute `./a25 -v`, `./s25 -v`, `./v25 -v` or `./525 -v`

NB: Omit the **-v** option if planning to use the `time` shell call

I've set the default compiler to `clang-13` in the Makefile, as this produces
the fastest executables in my testing, but edit the Makefile and switch to
`gcc` if that's what you have

a25, s25, v25, and 525 all support the **-v** option which will emit run-time metrics
If measuring with the `time` shell call, be aware that it typically takes about 1ms
for an executable to get up and running, the times reported by `time` will be slightly
higher than the internally measured times.  Also, the time taken to print the metrics
can become a non-neglible factor at these low execution times.

For speed, all solutions are written to a file named `solutions.txt` in the
current directory

`[a25|s25|v25|525] [-v] [-t num_threads] [-f word-file]`

- **-v** : Normally no console output is produced.  `-v` allows the executable to emit metrics
- **-t** : Allows the user to specify the number of threads to use.  By default the executables will use 1 or 2 less threads than there are CPUs on the system
- **-f** : Allows the user to specify an input word file to use.  By default the executables will use the words-alpha.txt file


### Execution Times

My development systems are a desktop AMD 5950x based PC, and an Intel i7-1165G7
based laptop. On my AMD system, which runs at a fixed at 4.7GHz CPU clock with
hyper-threading disabled, the following (internal) times are seen, which include
the loading of the full 4MB words_alpha.txt file:

**a25** takes around 18ms to complete, using 16 threads

**s25** takes around 0.98ms to complete using 16 threads

**v25** takes around 0.82ms to complete using 16 threads

**525** won't run on my AMD system.  On my Intel laptop it takes 1.35ms to complete
using 8 threads.  I estimate that on a full desktop runtimes of 0.80ms should be achievable

All algorithms use a bit-wise representation of the words for efficiency of comparing

All algorithms also implement lockless thread co-ordination for speed.


### History

My first attempt at a solution (source not included) focused on the New York Times
wordle list, for which there is 10 valid solutions.  It split the word set into
2-vowel and 1-vowel sub-sets, and utilised the fact that all actual words must have
one vowel in them (vowel being a, e, i, o , u, y).  This meant that there could be
either 5 words with 1-vowel, or 4 x 1-vowel words and a single 2-vowel word as
valid solutions.  That code just did a straight reduction search on the 1-vowel
set, which for the NYT Wordle List, only meant a little over 1-thousand 1-vowel
words.  The single 2-vowel set was then run over all 4 x 1-vowel sets to find
solutions.  This ran extremely quickly, being able to find all 10 solutions in
around 1.5ms

After commenting on the Youtube video, I then discovered that most people were
using the full words-alpha.txt file as a starting set, and even crazier to me,
people were using the full ~4MB word set, instead of at least using just a set
of only the 5 letter words. Okay, I guess?

My first solution did not work properly on the words-alpha input set as that set
contains many words with no-vowels in it, thereby breaking my assumption that
all valid words should have at least 1 vowel.

So I needed to redesign


### a25

This was my 2nd attempt at a full solution.  It implements a more generic version
of my first attempt.  Instead of using 2-vowel and 1-vowel words, it split the
input set into one set containing the most frequently character, which also
handily happens to be `a`, and everything else.  Using my first approach, this
generates a set of 4-word partial solutions from the non-a words.  For
completeness it also tries to find 5-word solutions using all letters but `a`,
but none exist.

Once all 4-word partial solutions are found, the `a` set is combined with the
4-word solutions to find all 5-word solutions.

This approach works very well and is very fast, but the idea of splitting the
input into multiple sets grew into my third solution attempt.


### s25

I've included `combos.c` which is non-functional pseudo-code that formed the
basis of `s25`.  I wrote a frequency implementation that originally ordered
all words from most the frequently occurring letter to the least.  This
approach produced just 9 word sets to represent all words and was easy to
combinate over.

After putting it all together, it worked, but quite slowly.  Taking about 50s
to run.  I then remembered some of the comments in the Youtube video about
trying from least-frequent to most-frequent, so I reversed the sort order, and
the solution now completed in about 15s.  This was still disappointing.

I then realised that we only need to use 25 letters out of 26, and so I tried
stopping the combination sequence 1 set short of the full combination set.
Doing this got the execution time down to under 1s while still finding all 538
solutions.

After this, I starting looking at other solutions online, and came across
Sylvester Hesp's solution, and to my surprise, we had both essentially developed
the same solution.  Where I was focusing on combinating over sets, Sylvester was
combinating over letters, but the actual implementation was extremely similar
in a remarkable example of convergent design.

Sylvester's work is here: https://github.com/oisyn/parkerwords

I then saw how Sylvester had solved the 25 instead of 26 approach with his
clever skipping implementation that skips at the start of the combinator,
instead of my approach that skipped at the end.  His approach actually
applied perfectly to my code (I owe a huge debt of gratitude to his solution
to this issue), and with that the initial implementation of `s25` was running
in around 4-5ms for the actual main algorithm loop, and ~5-7ms to load in the
the words-alpha file, with overall run-times ranging from 9-14ms depending on
what the system is doing.  Run times were signiciantly less if using the
words-alpha-five.txt or NYT files as input.


### v25

`v25.c` is just `s25.c` but with AVX2 instructions added for processing the
key sets.  The use of AVX-2 saw the main algorithm loop be sped up by almost
3x.  This is my first time ever playing around with AVX instructions, so it
was a fun learning experience.


### 525

`525.c` is `v25.c` with AVX-512 instructions, instead of AVX2.  This will not
compile successfully on any platform that does not support AVX-512.

AVX-512 will run the main algorithm in around 65% of the time of the AVX2
version.  Of course, file load times are still a significant factor, and so
I'd expect to see overall run-times at ~75-80% of `v25` run times.  Sadly I
do not have a high-end AVX-512 capable desktop to verify this claim, although
these ratios were true on the AVX-512 laptop.


### Words Alpha File Reading

A large problem I dealt with was how to quickly read in and process the 4MB
word file.  My first attempt just processed the file with a straight mmap()
and that saw the file get read in and processed in around 7ms.  My second
attempt used multiple threads, which each thread working on their own portion
of the input file.  This worked fairly well, but hash table building was a
large source of contention. Additionally, the system doesn't really like
many threads jumping about all over a single file.

My final solution was to set the readers to incrementing an atomic integer to
gain access to the next 8KB block of the file.  This improved locality and
also ensured a (mostly) sequential step through the file, which boosted file
reading times dramatically.  I then had the readers focus on finding 5 letter
words and inserting them into an array, with an atomic integer defining where
in the array any particular reader would insert their word.

I then had the main thread process the word array that was being built by the
reader threads to build the hash table concurrently.  Doing all this together
managed to get file load and hash table build times to under 0.7ms on my AMD
system.


### Frequency Rescanning

About mid-way through the frequency set build, we rescan the frequencies and
re-order the remaining sets from most-frequent to least.  The works because
the bulk of the benefit from the least to most frequent ordering has already
been gained by this stage, and reducing the length of the combinatorial "tail"
has an amortizing effect of more than making up for the increased set sizes
mid-algorithm


### Clean-Up

I moved all the file reading and other utility functions to a common `utilities.h`
file that is included by the main algorithms.  In this way, any future tweaks
should be consistent across all implementations

### Major Update Sep 1st

I received an email from Landon Kryger about his solution here:
https://github.com/GuiltyBystander/5words

Landon had come up with the most efficient implementation of the search problem
I'd ever seen.  Truly fantastic work!  Landon says that he arrived at his
solution also fairly independently, but at its core it uses the same basic
algorithm that Sylvester and myself had found.

Where Landon went one step further is that he had added a further breakdown of
the search space, where the presence (or not) of the most common letters is
used to divide each set further, such that when progressing through the sets,
characters that are already seen are used to eliminate sub-sets of the sets
from being considered.  This comes at a cost of a somewhat high amount of set
up overhead, but for a single-threaded solution, his was simply the fastest
around.

Landon was looking for someone to collaborate with to help make his solution
faster.  I spent some time analysing how his algorithm works, and while it is
extremely clever and elegant, it is somewhat difficult to parallelise the
setup steps, and the bit-mask remapping he implements which is essential to
the speed of his approach consumes a significant amount of CPU time, and also
forces a point of linearity to derive the mapping.  The remapping of the bitmaps
can be parallelised, as can the set building, but spatial locality issues arise
from the data being scattered across many different arrays/vectors.

As a direct comparison, even though Landon's solution makes 30x less comparisons
during the search phase, his algorithm is just 3x faster when using a single
thread.  I think we can implement a solution true to his original vision in a
highly parallelised manner, but it would be a fair amount of work.

So instead I worked with Landon to arrive at a hybrid solution that isn't as
elegant or as flexible as Landon's core solution, but preserves the strong
spatial locality inherent in the "pointers within a single set" approach that
my solution implements.

With that, I updated s25.c, v25.c and 525.c to implement this hybrid approach,
and huge speedups were seen in low thread numbers, but smaller speedups as
threads are increased.

v25.c is now capable of finding all 538 solutions on my system in under 1ms
including loading the 4MB words-alpha.txt file, solving, and writing the
results out.  Even using just a single thread takes under 6ms.

Where Landon's approach uses the 6 most frequently occurring characters to
create subsets from, the hybrid solution only uses 2, and it is still slightly
slower than Landon's solution for a single thread scenario for non-AVX mode.

Update:  The `master` branch has 6 tier code now

I was able to parallelise the tier setup paths.  This has resulted in the setup
for 6 tiers taking only 10us more than the setup for 3 tiers.  As a result, the
6 tier path is now the fastest in all scenarios.

Update 2:  I discovered a low-CPU-cost way to calculate the best characteres to
use to split subsets for each character set.  Main algorithm compares were
reduced by 33%.


### Conclusion

I'm including the words-alpha file, plus a version of words-alpha that only has
the 5-letter words in it, as well as the NYT Wordle sets.  It seems silly to me
for us to be spending time finding ~80KB of 5 letter words from a ~4MB file but
since that's what people seem to be doing, so I included the full set here.

Thank you to Matt Parker for making the problem public, and all the intelligent
and robust discussion of solutions in the Youtube comments, and a big thank you
to Sylvester Hesp for his approach to the set skipping issue, and to Landon
Kryger for his insightful algorithm improvements!

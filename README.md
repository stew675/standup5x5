### Stand Up Maths 5 x 5 unique letters solutions

My solutions to the problem posed by Matt Parker of Stand Up Maths

See: https://www.youtube.com/watch?v=_-AfhLQfb6w


### Building and Running

Just run make and execute `./a25` or `./s25`

**a25** should take around 0.025-0.050s to complete depending on your hardware

**s25** should take around 0.008-0.016s to complete depending on your hardware

Both algorithms uses a bit-wise representation of the words for efficiency of comparing

`s25` is the more recent solution, and takes some arguments.  For speed, all
solutions are written to a `solutions.txt` file.

The arguments are:

`s25 [-v] [-t num_threads] [-f word-file]`

- _-v_ : `s25` normally produces no console output.  `-v` allows it to emit metrics
- _-t_ : Allows the user to specify the number of threads to use.  By default `s25` will used 2 less threads than there are CPUs on the system
- _-f_ : Allows the user to specify an input word file to use.  By default `s25` will use the alpha words set.


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
of only 5 letter words. Okay, I guess?

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
completes it also tries to find 5-word solutions using all letters but `a`,
but none exist.

Once all 4-word partial solutions are found, the `a` set is combined with the
4-word solutions to find all 5-word solutions.

This approach works very well and is very fast, but the idea of splitting the
input into multiple sets grew into my third solution attempt.


### s25

You can see some dead combination code in `a25.c` which was used to create the
idea that we could combinate over a small number of word sets to derive all
solutions.  I included `combos.c` which is non-functional pseudo-code that
formed the basis of `s25.  I wrote a frequency implementation that originally
ordered all words from most frequently occurring letter to the least.  This
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

I then saw how Sylvester had solved the 25 instead of 26 approach with his
clever skipping implementation that skips at the start of the combinator,
instead of my approach that skipped at the end.  His approach actually
applied perfectly to my code (I owe a huge debt of gratitude to his solution
to this issue), and with that `s25` was running in around 4-5ms for the actual
main algorithm loop, and ~5-7ms to load in the the words-alpha file, with
overall run-times ranging from 9-14ms depending on what the system is doing.
Run times are much less if using the words-alpha-five.txt or NYT files as input.

Both a25.c and s25.c implement lockless thread co-ordination for speed.


### v25

`v25.c` is just `s25.c` with AVX2 instructions add for processing the key sets.
Where `s25.c` takes ~4ms to process the main algorithm loop, `v25.c` will does
it in around 2.1ms on my system.  This is my first time ever playing around
with AVX instructions, so it was a fun learning experience.


### 525

`525.c` is `v25.c` with AVX-512 instructions, instead of AVX2.  This will not
compile successfully on any platform that does not support AVX-512.

AVX-512 will run the main algorithm in around 2/3rds of the time of the AVX2
version.  Of course, file load times are still a significant factor, and so
overall expect to see overall run-times at ~75% of `v25` run times


### Conclusion

I'm including the words-alpha file, plus a version of words-alpha that only has
the 5-letter words in it, as well as the NYT Wordle sets.  It seems silly to me
for us to be spending time finding ~80KB of 5 letter words from a ~4MB file but
since that's what people seem to be doing, I'll include the full set here as well

Thank you to Matt Parker for making the problem public, and all the intelligent
and robust discussion of solutions in the Youtube comments, and a big thank you
to Sylvester Hesp for his approach to the set skipping issue.

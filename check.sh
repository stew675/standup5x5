echo
echo "SOLUTION VALIDATOR"
echo
echo "Note that some differences may occur as this is dependent upon"
echo "which permutation of a word makes it into the hash-table first"
echo "This can occur with many concurrent reader threads"
echo
echo "Inspect the diff output to ensure that the only differences"
echo "are those with just the characters re-arranged"

echo
echo
echo "Checking s25 output correctness"
rm -f solutions.txt
./s25 -f words_alpha.txt
sort < solutions.txt | diff - expected_solutions.txt

echo
echo
echo "Checking v25 output correctness"
rm -f solutions.txt
./v25 -f words_alpha.txt
sort < solutions.txt | diff - expected_solutions.txt

echo
echo
echo "Checking 525 output correctness"
rm -f solutions.txt
./525 -f words_alpha.txt
if [ -e "solutions.txt" ]; then
	sort < solutions.txt | diff - expected_solutions.txt
fi


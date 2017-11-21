#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/basic.h>
#include <iostream>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <vector>
#include <string>
#include <map>

#include "w50_functions.h"

using namespace seqan;
using namespace std;

// Does the rounding work correctly
SEQAN_DEFINE_TEST(rounding_test_1)
{
   SEQAN_ASSERT_EQ(roundUp(25,50), 50);
}

// Does the rounding work when it's directly on the top edge of the bin
SEQAN_DEFINE_TEST(rounding_test_2)
{
   SEQAN_ASSERT_EQ(roundUp(1000,50), 1000);
}

// Does the rounding work when it's directly at the bottom
SEQAN_DEFINE_TEST(rounding_test_3)
{
   SEQAN_ASSERT_EQ(roundUp(950,50), 950);
}

// Does the check function return 1 (failure) when out of order 
// reference is encountered
SEQAN_DEFINE_TEST(check_sorted_test_1)
{
   vector<CharString> chromosomes;
   chromosomes.push_back("chr1");
   chromosomes.push_back("chr2");
   chromosomes.push_back("chr3");
   CharString current = "chr2";

   // Only required to print out the file names
   ModifyStringOptions options;

   SEQAN_ASSERT_EQ(checkSorted(chromosomes, current, options, false), 1);
}

// Does the check function return 0 (success) when correct order
SEQAN_DEFINE_TEST(check_sorted_test_2)
{
   vector<CharString> chromosomes;
   chromosomes.push_back("chr1");
   chromosomes.push_back("chr2");
   chromosomes.push_back("chr3");
   CharString current = "chr4";

   // Only required to print out the file names
   ModifyStringOptions options;

   SEQAN_ASSERT_EQ(checkSorted(chromosomes, current, options, false), 0);
}

// 
SEQAN_DEFINE_TEST(check_insertIntoMap_1)
{

//insertIntoMap(std::map<int, WindowValues> &bins,
  //                ModifyStringOptions options, int &largest,
    //              CharString &currentRef, GffRecord &record)
}

SEQAN_BEGIN_TESTSUITE(test_w50_creator)
{
   // Call roundUp() tests
   SEQAN_CALL_TEST(rounding_test_1);
   SEQAN_CALL_TEST(rounding_test_2);
   SEQAN_CALL_TEST(rounding_test_3);

   // Call checkSorted() tests
   SEQAN_CALL_TEST(check_sorted_test_1);
   SEQAN_CALL_TEST(check_sorted_test_2);
}
SEQAN_END_TESTSUITE

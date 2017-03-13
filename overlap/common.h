#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/stream.h>    // to stream a CharString into cout
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <math.h>       /* sqrt */
#include <seqan/store.h> /* FragmentStore */
#include <queue>
#include <vector>
#include <ctime>
#include "boost/multi_array.hpp"
#include <cassert>
#include <boost/unordered_map.hpp>
#include <string>
#include <thread>
#include <mutex>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <iomanip>
#include <functional>
#include "boost/icl/interval.hpp"
#include "boost/icl/interval_map.hpp"
#include <boost/icl/split_interval_map.hpp>
#include <boost/icl/interval_set.hpp>
#include <boost/algorithm/string.hpp>
#include "boost/lexical_cast.hpp"
#include <set>

#include <chrono>
using namespace seqan;
using namespace std;

using namespace boost::icl;
using boost::bad_lexical_cast;

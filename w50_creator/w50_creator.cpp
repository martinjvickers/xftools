#include "common.h"

struct ModifyStringOptions
{
        CharString inputFileName;
	int window_length;
};

struct WindowValues
{
        int c;
        int t;
	int n;
};

seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
	seqan::ArgumentParser parser("w50_creator");
	addOption(parser, seqan::ArgParseOption("i", "input-file", "Path to the input file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "input-file");
	setShortDescription(parser, "Methylation Tools");
	setVersion(parser, "0.0.1");
	setDate(parser, "July 2016");
	addUsageLine(parser, "-i sequence.fastq [\\fIOPTIONS\\fP] ");
	addOption(parser, seqan::ArgParseOption("l", "window-length", "Size of window",seqan::ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "window-length", "50");

	addDescription(parser, "Create a w50 file from a w1 file.");
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res;

	getOptionValue(options.inputFileName, parser, "input-file");
	getOptionValue(options.window_length, parser, "window-length");

	return seqan::ArgumentParser::PARSE_OK;

}

int roundUp(int numToRound, int multiple)
{
	if (multiple == 0)
		return numToRound;

	int remainder = numToRound % multiple;
	if (remainder == 0)
		return numToRound;

	return numToRound + multiple - remainder;
}

/*
Aim: Safestyle, BOGOF: using hashes to window a GFF

Current progress: It compiles!

*/
int main(int argc, char const ** argv)
{

	//parse our options
	ModifyStringOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
	
	GffFileIn gffFileIn;
	if (!open(gffFileIn, toCString(options.inputFileName)))
	{
		std::cerr << "ERROR: Could not open the file.\n";
		return 1;
	}

	CharString currentRef = NULL;

	//define our map to store windows
	//typedef boost::unordered_map<int, WindowValues> map;
	map<int, WindowValues> map;

	std::cout << std::setprecision(4) << std::fixed;

	//stores chromosomes we've seen
	vector<CharString> haveSeen;

	int count = 0;

	// Copy the file record by record.
	GffRecord record;
	while (!atEnd(gffFileIn))
	{

		try
		{
			readRecord(record, gffFileIn);
			count++;
		}
		catch (Exception const & e)
		{
			std::cerr << "ERROR: " << e.what() << std::endl;
			std::cerr << "ERROR: " << record.ref << " " << currentRef << " "<< count << std::endl;
			return 1;
		}

		//First Chromosome
		if(currentRef == NULL)
		{
			currentRef = record.ref;		
		}

		//at this point we will have dump the contents of the hash
		if(record.ref != currentRef)
		{
			
			//search to see if we've already seen this chromosome
			for(auto const& value: haveSeen) 
			{
				if(value == currentRef)
				{
					cerr << "You're about to write out " << currentRef << " but we've already seen this chromosome. This probably means your w1 file is not in chromosome order."<< endl;
					cerr << "To resolve this, simply sort your file by Chromosome order, something like; sort -kn 1 " << toCString(options.inputFileName) << " > " << toCString(options.inputFileName) << "_sorted.w1.gff"<< endl;
				}
			}

			for(auto& p: map)
			{
				float score;
				//write out GFF
				if ((p.second.c > 0) || (p.second.t > 0)){
					score = (float)p.second.c / (float)(p.second.c + p.second.t);
				} else {
					score = 0.0000;
				}
					cout << currentRef << "\t.\twindow\t"<< p.first - (options.window_length-1) << "\t" << p.first << "\t"<< score << "\t.\t.\tc=" << p.second.c << ";t=" << p.second.t << ";n=" << p.second.n << endl;
			}
	
			//we should save the current reference as a check to see if we come across it again out of order
			haveSeen.push_back(currentRef);

			//clear map
			map.clear();

			currentRef = record.ref;
		}

		//here, we add to the hash
		if(record.ref == currentRef)
		{

			//////ARGH!!!
			/*
			SeqAn uses 0-based half-open intervals for internal representation
			however all the previous dzlab stuff (at least w50 creation so far)
			uses a 1-based half-open interval.

			READ HERE
			https://seqan.readthedocs.io/en/master/Tutorial/InputOutput/GffAndGtfIO.html
			*/
			int window = roundUp(record.beginPos+1, options.window_length);

			//we expect the following format of tags and values  c=4;t=0;n=1
			//but we will not assume they are always in that order when reading 
			//them from the file. Since we have a struct to store the values,
			//we don't need to store the tag names themselves, just values.
			
			int c;
			int t;
			int n;
			
			WindowValues currWindowVal;

			bool haveN = 1;

			//iterate through looking for c,t,n (will ignore anything else)
			for(int i = 0; i < length(record.tagNames); i++)
			{	
				
				if(record.tagNames[i] == "c")
				{
					currWindowVal.c = atoi(toCString(record.tagValues[i]));
				}

                                if(record.tagNames[i] == "t")
                                {
                                        currWindowVal.t = atoi(toCString(record.tagValues[i]));
                                }

                                if(record.tagNames[i] == "n")
                                {
                                        currWindowVal.n = atoi(toCString(record.tagValues[i]));
					haveN = 0;
                                }

				//need to check if we have all three??

			}

			//means we have no n's in the input file, so we add one
			if(haveN == 1)
			{
				currWindowVal.n = 1;
			}

			//get element in hash table, if not there, just put this one in it
			//if window exists, get values and add current Window to it
			WindowValues get = map[window];
			currWindowVal.c = get.c + currWindowVal.c;
			currWindowVal.t = get.t + currWindowVal.t;
			currWindowVal.n = get.n + currWindowVal.n;

			//now insert back into map
			map[window] = currWindowVal;

			//update currentRef
			currentRef = record.ref;

		}
	}

	//read out very last chromosome
	for(auto& p: map)
	{

		//search to see if we've already seen this chromosome
                for(auto const& value: haveSeen)
                {
                	if(value == currentRef)
                        {
                        	cerr << "You're about to write out " << currentRef << " but we've already seen this chromosome and written it out to file. This probably means your w1 file is not sorted in chromosome order."<< endl;
                                cerr << "To resolve this, simply sort your file by Chromosome order, something like; sort -kn 1 " << toCString(options.inputFileName) << " > " << toCString(options.inputFileName) << "_sorted.w1.gff"<< endl;
                        }
                }


		float score;
		//write out GFF
		if ((p.second.c > 0) || (p.second.t > 0)){
			score = (float)p.second.c / (float)(p.second.c + p.second.t);
		} else {
			score = 0;
		}
		cout << currentRef << "\t.\twindow\t"<< p.first - (options.window_length-1) << "\t" << p.first << "\t"<< score << "\t.\t.\tc=" << p.second.c << ";t=" << p.second.t << ";n=" << p.second.n << endl;
	}


	return 0;
}

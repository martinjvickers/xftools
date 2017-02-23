#include "common.h"
#include <unordered_map>
#include "key.cpp"
#include <functional>
#include "boost/icl/interval.hpp"
#include "boost/icl/interval_map.hpp"
#include <boost/icl/split_interval_map.hpp>
#include <boost/icl/interval_set.hpp>
#include <set>
using namespace boost::icl;

struct ModifyStringOptions
{
        CharString inputFileName;
	CharString inputAnnotationFileName;
	CharString outputFileName;
	bool exclude;
};

class Feature
{
public:
	Feature():_score(0),_strand('*'),_ref("Meh"),_endPos(0),_startPos(0){} // default constructor which i hope to never use?
	Feature(int score, char strand, CharString ref, int startPos, int endPos):_score(score),_strand(strand),_ref(ref),_endPos(endPos),_startPos(startPos){} // constructor
	int endPos()const  {return _endPos;}
	int startPos()const  {return _startPos;}
	int score()const{return _score;}
	CharString ref()const{return _ref;}
	char strand()const{return _strand;}

	void increment(){_score++;}

	//next to implement
		//give the tags(e.g. c,t,n)

	Feature& operator += (const Feature& right)
	{	
		_score += right.score(); 
		return *this; 
	}

private:
	int _startPos;
	int _endPos;
	char _strand;
	CharString _ref;
	int _score;
};

bool operator == (const Feature& left, const Feature& right)
{
	//== if everything (except for score) is the same 
	return left.strand()==right.strand() && left.ref()==right.ref() && left.endPos()==right.endPos() && left.startPos()==right.startPos(); 
}

//map of map of split interval map. 
typedef split_interval_map<int, Feature> featuremap;
typedef std::map<char, featuremap> strandedmap;
typedef std::map<string, strandedmap> completemap;

seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
	seqan::ArgumentParser parser("wba");
	addOption(parser, seqan::ArgParseOption("i", "input-file", "Path to the input file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "input-file");
	addOption(parser, seqan::ArgParseOption("a", "input-annotation-file", "Path to the input filter file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
        setRequired(parser, "input-annotation-file");
	addOption(parser, seqan::ArgParseOption("o", "output-file", "Path to the output file", seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
	setRequired(parser, "output-file");

	setShortDescription(parser, "Window by Annotation");
	setVersion(parser, "0.0.1");
	setDate(parser, "Feb 2017");
	addUsageLine(parser, "-i input.gff -a annotation.gff -o output.gff [\\fIOPTIONS\\fP] ");

	addDescription(parser, "Given an annotation/feature file and an input w1 file, give counts for each feature");
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res;

	getOptionValue(options.inputFileName, parser, "input-file");
	getOptionValue(options.inputAnnotationFileName, parser, "input-annotation-file");
	getOptionValue(options.outputFileName, parser, "output-file");

	return seqan::ArgumentParser::PARSE_OK;
}

/*
Aim: Given an annotation file and an input w1 file, this will add all of the w1 records that are within that annotation file.
Current progress: It works in principle
*/
int main(int argc, char const ** argv)
{
	//parse our options
	ModifyStringOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	split_interval_map<int, Feature> map; // old way
	completemap results; // new way results[ref][strand].thingy

	//get our value and add value to map
	GffFileIn gffAnnotationIn;
	GffRecord annotationrecord;
        if (!open(gffAnnotationIn, toCString(options.inputAnnotationFileName)))
        {
                std::cerr << "ERROR: Could not open example.gff" << std::endl;
                return 1;
        }

	//put all of our annotations into a vector and store the pairs in an interval map
	while (!atEnd(gffAnnotationIn))
        {
		readRecord(annotationrecord, gffAnnotationIn);
		Feature feat(0, annotationrecord.strand, annotationrecord.ref, annotationrecord.beginPos, annotationrecord.endPos);
		discrete_interval<int> inter_val = discrete_interval<int> (annotationrecord.beginPos, annotationrecord.endPos);
		map += make_pair(inter_val, feat);
		//let's allso append to the correct place in the new completemap
		results[toCString(annotationrecord.ref)][annotationrecord.strand] += make_pair(inter_val, feat);
	}

	close(gffAnnotationIn);

	//now load our new data
	GffFileIn gffRawIn;
	GffRecord to_bin_record;
	if (!open(gffRawIn, toCString(options.inputFileName)))
	{
		std::cerr << "ERROR: Could not open example.gff" << std::endl;
                return 1;
	}

	//logic, this is where we bin or count etc with our overlaps.
	while (!atEnd(gffRawIn))
        {
		readRecord(to_bin_record, gffRawIn);

		discrete_interval<int> key = discrete_interval<int> (to_bin_record.beginPos, to_bin_record.endPos);
		std::pair<split_interval_map<int, Feature>::iterator, split_interval_map<int, Feature>::iterator> itres = map.equal_range(key);

		//this is the old way
		for(auto it = itres.first; it != itres.second; ++it)
			(*it).second.increment();

		//without wrecking the old way, let's increment the new complete map


		//check if stranded
		std::pair<split_interval_map<int, Feature>::iterator, split_interval_map<int, Feature>::iterator> itresnew;
		if(results[toCString(to_bin_record.ref)][to_bin_record.strand].size() > 0)
		{
                	itresnew = results[toCString(to_bin_record.ref)][to_bin_record.strand].equal_range(key);
			for(auto itnew = itresnew.first; itnew != itresnew.second; ++itnew)
	                        (*itnew).second.increment();
		} else {
			//cout << "Input data has no strand information, so adding to both strands" << endl; //probably only show this message once
			for(auto& i : results[toCString(to_bin_record.ref)])
			{
				itresnew = i.second.equal_range(key);
				for(auto itnew = itresnew.first; itnew != itresnew.second; ++itnew)
					(*itnew).second.increment();
			}
		}
		
		//DO NOT DELETE. THIS IS THE BASE LOGIC OF THIS.
//		std::pair<split_interval_map<int, Feature>::iterator, split_interval_map<int, Feature>::iterator> itresnew = results[toCString(to_bin_record.ref)][to_bin_record.strand].equal_range(key);
//		for(auto itnew = itresnew.first; itnew != itresnew.second; ++itnew){
//                        (*itnew).second.increment();
//		}
	}

	close(gffRawIn);

	//now, let's reopen out feature file and get our feature counts
	GffFileIn gffFeatureIn;
	GffRecord featurerecord;
	if (!open(gffFeatureIn, toCString(options.inputAnnotationFileName)))
        {
                std::cerr << "ERROR: Could not open example.gff" << std::endl;
                return 1;
        }

	while (!atEnd(gffFeatureIn))
        {
		readRecord(featurerecord, gffFeatureIn);
		discrete_interval<int> key = discrete_interval<int> (featurerecord.beginPos, featurerecord.endPos);
		std::pair<split_interval_map<int, Feature>::iterator, split_interval_map<int, Feature>::iterator> itres = map.equal_range(key);

		int score = 0;
		for(auto it = itres.first; it != itres.second; ++it)
                {
			Feature mmm = (*it).second;
			score = score + mmm.score();
                }


		//again, without wrecking the old way, lets not extract our our features from the completemap and print them
		score = 0;
		std::pair<split_interval_map<int, Feature>::iterator, split_interval_map<int, Feature>::iterator> itresnew = results[toCString(featurerecord.ref)][featurerecord.strand].equal_range(key);
		for(auto itnew = itresnew.first; itnew != itresnew.second; ++itnew)
                {
			Feature mmm = (*itnew).second;
			score = score + mmm.score();
		}
		cout << featurerecord.ref << "\t" << featurerecord.source << "\t" << featurerecord.type << "\t" << featurerecord.beginPos << "\t" << featurerecord.endPos << "\t" << score << "\t" << featurerecord.strand << "\t" << featurerecord.phase << endl;
	}

	close(gffFeatureIn);

	return 0;
}

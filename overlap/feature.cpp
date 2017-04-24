
class Feature
{
public:
        Feature():_score(0),_strand('*'),_ref("Meh"),_endPos(0),_startPos(0){} // default constructor which i hope to never use?
        Feature(int score, char strand, CharString ref, int startPos, int endPos):_score(score),_strand(strand),_ref(ref),_endPos(endPos),_startPos(startPos){} // constructor
        int endPos()const  {return _endPos;}
        int startPos()const  {return _startPos;}
        int score()const{return _score;}
        std::map<string, string> tagmap()const{ return _tagmap;}
        CharString ref()const{return _ref;}
        char strand()const{return _strand;}

        void increment()
        {
                _score++;
        }

        //the problem here is that you can put anything into a tag. So, we will do all the logic in this function
        //a methyl tag looks like c=0;t=4
        // i want a generic method for this. as in, if it's an integer or float, add it, if it's a string, just ignore.
        //not sure how to do this though. Let's start by simply doing what we know. c's and bloody t's.
        ///////AND WHY ARE STRING SETS SETS? WHY NOT MAPS?!
        void addtags(StringSet<CharString> tagNames, StringSet<CharString> tagValues)
        {
                for(int i = 0; i < length(tagNames); i++)
                {
                        //for each tag, find it in the _tagmap
                        std::map<string, string>::iterator it;
                        it = _tagmap.find(toCString(tagNames[i]));
                        //if it exists, do something
                        if(it != _tagmap.end())
                        {
                                double toadd, current;
                                current = 0;
                                string currentstr = _tagmap[toCString(tagNames[i])];
                                try {
                                        toadd = boost::lexical_cast<double>(toCString(tagValues[i]));
                                        current = boost::lexical_cast<double>(_tagmap[toCString(tagNames[i])]);
                                }
                                catch (bad_lexical_cast &) {
                                        //don't actually do anything then
                                }
                                _tagmap[toCString(tagNames[i])] = std::to_string(toadd + current);
                        }
                        else //just add it
                        {
                                _tagmap[toCString(tagNames[i])] = toCString(tagValues[i]);
                        }
                }
        }

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
        std::map<string, string> _tagmap;
};

bool operator == (const Feature& left, const Feature& right)
{
        //== if everything (except for score) is the same 
        return left.strand()==right.strand() && left.ref()==right.ref() && left.endPos()==right.endPos() && left.startPos()==right.startPos();
}

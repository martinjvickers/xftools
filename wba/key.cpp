class MyPair
{
	int b_value;
	int e_value;	
	CharString id;

	public:
		MyPair(CharString refid, int begin, int end)
		{
			id = refid;
			b_value = begin;
			e_value = end;
		};
		friend bool operator< (const MyPair &c1, const MyPair &c2);
		friend bool operator> (const MyPair &c1, const MyPair &c2);
		friend bool operator== (const MyPair &c1, const MyPair &c2);

		int getBegin() const
		{
		        return b_value;
		};
		int getEnd() const
		{
		        return e_value;
		};
		CharString getID() const
                {
                        return id;
                };
		bool within(CharString val_id, int val) const
		{
			string c1str = toCString(id);
			string c2str = toCString(val_id);
			cout << c1str << " " << c2str <<" " << c1str.compare(c2str) << endl;
			if(c1str.compare(c2str) && val >= b_value && val <= e_value)
				return true;
			else
				return false;
		};
};

bool operator< (const MyPair &c1, const MyPair &c2)
{
        string c1str = toCString(c1.id);
        string c2str = toCString(c2.id);
	//cout << c1str << " " << c2str <<" " << c1str.compare(c2str) << endl;
        if(c1str.compare(c2str) && c1.e_value < c2.e_value && c1.b_value < c2.e_value)
                return true;
        else
                return false;
}

bool operator> (const MyPair &c1, const MyPair &c2)
{
        string c1str = toCString(c1.id);
        string c2str = toCString(c2.id);
	//cout << c1str << " " << c2str <<" " << c1str.compare(c2str) << endl;
	if(c1str.compare(c2str) && c1.e_value > c2.e_value && c1.b_value > c2.e_value)
		return true;
	else
		return false;
}

bool operator== (const MyPair &c1, const MyPair &c2)
{
	string c1str = toCString(c1.id);
	string c2str = toCString(c2.id);
	//cout << c1str << " " << c2str <<" " << c1str.compare(c2str) << endl;
        if(c1str.compare(c2str) && c1.e_value == c2.e_value && c1.b_value == c2.b_value)
                return true;
        else
                return false;
}



/*
bool operator<= (const MyPair &c1, const MyPair &c2)
{
        return c1.e_value <= c2.b_value;
}

bool operator>= (const MyPair &c1, const MyPair &c2)
{
        return c1.e_value >= c2.b_value;
}
*/

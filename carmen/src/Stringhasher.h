#ifndef HASH_TYPEDEFS
typedef boost::unordered_map<std::string, int> HASH_S_I;
typedef boost::unordered_map<std::string, std::string> HASH_S_S;
#define HASH_TYPEDEFS 1
#endif

class stringhasher{
public:
	static const int EQ=0;
	static const int LT=-1;
	static const int GT=1;

	static void stringhasher::split(std::string line, std::vector<std::string>* arr, char* token){
		arr->clear();
		typedef boost::tokenizer< boost::char_separator<char> > tokenizer;
		boost::char_separator<char> sep(token);
		tokenizer tok(line, sep);
		for(tokenizer::iterator iter=tok.begin(); iter!=tok.end(); iter++){
			arr->push_back( (*iter) );
		}
	}

	static void stringhasher::join(std::string& result, std::vector<std::string>* arr, std::string token){
		result = "";
		int N =(int)arr->size();
		for(int i=0; i<N; i++){
			result += arr->at(i);
			if(i<N-1)
				result += token;
		}
	}

	static void stringhasher::get_idval(std::string rulepart, std::string& identifier, std::string& value){
		int loc_is = (int)rulepart.find("_is_");
		int len = (int)rulepart.length();
		if(loc_is>-1){
			identifier = rulepart.substr(0,loc_is);
			value = rulepart.substr(loc_is+4, len-loc_is);
		}
		else
			throw new std::string("Rule part does not contain _is_, cannot parse it.");
	}

	static void stringhasher::get_id_dir_val(std::string str, std::string& identifier, int& dir, float& value){
		int loc_LT = (int)str.find("<");
		int loc_GT = (int)str.find(">");
		int loc_EQ = (int)str.find("=");
		int len = (int)str.length();
		if( loc_EQ != -1 ){
			dir = stringhasher::EQ;
			identifier = str.substr(0,loc_EQ);
			value = boost::lexical_cast<float>( str.substr(loc_EQ+1, len-loc_EQ ) );
		}
		else if( loc_LT != -1 ){
			dir = stringhasher::LT;
			identifier = str.substr(0,loc_LT);
			value = boost::lexical_cast<float>( str.substr(loc_LT+1, len-loc_LT ) );
		}
		else if( loc_GT != -1 ){
			dir = stringhasher::GT;
			identifier = str.substr(0,loc_GT);
			value = boost::lexical_cast<float>( str.substr(loc_GT+1, len-loc_GT ) );
		}
		else
			throw new std::string("Rule part does not contain <, >, or =, cannot parse it.");
	}

	static void stringhasher::parse_limit(std::string limit, std::string& str_attrib, std::string& str_value, bool& is_equal){
		// limit should have the form foo=bar. set str_attrib to foo and str_value to bar.
		int len = (int)limit.length();
		int loc_bang = (int)limit.find("!");
		if(loc_bang>-1){
			is_equal=false;
			str_attrib = limit.substr(0,loc_bang);
			str_value = limit.substr(loc_bang+1, len-loc_bang);
		}
		else{
			is_equal=true;
			int loc_equal = (int)limit.find("=");
			if(loc_equal==-1)
				throw std::string("Incorrect format for class limit: must include = or ! symbol.");
			str_attrib = limit.substr(0,loc_equal);
			str_value = limit.substr(loc_equal+1, len-loc_equal);
		}
	}

	size_t operator() (const std::string& s) const{
		size_t h = 0;
		std::string::const_iterator p, p_end;
		for(p = s.begin(), p_end = s.end(); p != p_end; ++p){
			h = 31 * h + (*p);
	    }
		return h;
	}
	bool operator() (const std::string& s1, const std::string& s2) const{
		return s1 < s2;
	}
};

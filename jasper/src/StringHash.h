/*
class stringhasher : public stdext::hash_compare <std::string>{
public:
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

typedef stdext::hash_map<std::string, int, stringhasher> HASH_S_I;
typedef stdext::hash_map<std::string, std::string, stringhasher> HASH_S_S;
*/
const bool OPT_REQUIRED = true;
const bool OPT_OPTIONAL = false;
class Option{
public:
	Option(std::string long_name, std::string short_name, std::string description, std::string default_value, bool required);
	std::string long_name;
	std::string short_name;
	std::string description;
	std::string value;
	std::string default_value;
	bool required;
	bool is_set;
	void set_value(std::string value);
	std::string get_value();
};

int read_args(int argc, char* argv[], std::vector<Option*>* options, std::string preamble);

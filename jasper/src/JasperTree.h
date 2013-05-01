enum{
	LEAF_ERROR = -1,
    LEAF_ROOT,
	LEAF_FILE,
	LEAF_QUERY,
	LEAF_SUBHEAD
};

class TreeData : public wxTreeItemData{
public:
	TreeData(std::string str, std::string original_label, int leaf_type);
	std::string str;
	std::string original_label;
	int leaf_type;
};


class JasperTree : public wxTreeCtrl{
public:
	JasperTree(wxWindow* parent, wxWindowID id, Investigation* investigation, wxSize size);
	void refresh_from_directory();
	void get_selected_label(std::string& label);
	void get_selected_original_label(std::string& label);
	void get_selected_filename(std::string& str);
	int get_selected_leaftype();
    void open_tree_to_file(std::string filename);

private:
	ConfigParser* cp;
	Investigation* investigation;
	wxWindow* parent;
	wxTreeItemId root;
	wxImageList imagelist;
	int img_folder;
	int img_file_open;
	int img_file_normal;

	void get_limit_and_label(std::string filename, std::string extension, std::string& limit, std::string& label);
	wxTreeItemId find_child_named(wxTreeItemId item_top, std::string query);
	wxTreeItemId add_label(wxTreeItemId item_q, std::string filename, std::string extension, std::string label);
	wxTreeItemId create_query(std::string label_to_show, std::string original_label);
	wxTreeItemId append_child_with_icons(wxTreeItemId parent, std::string label, TreeData* treedata);
};
